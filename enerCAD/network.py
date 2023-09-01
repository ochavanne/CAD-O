# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 18:14:49 2023

@author: Olivier Chavanne
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import os
from scipy.spatial import distance_matrix
from shapely.geometry import LineString, Point
import matplotlib.pyplot as plt
import datetime



def get_trace(zone_cad, gdf_centrale, power_EGID):
    
    df_DN = pd.read_csv("DN.csv", delimiter=";")
    df_sizing = df_DN.copy()
    gdf_cad = zone_cad.copy()

    # Perform the merge based on matching values between "RegBL_EGID" (gdf_cad) and "EGID" (Power_EGID)
    gdf_cad['RegBL_EGID'] = gdf_cad['RegBL_EGID'].astype('int64')
    power_EGID['EGID'] = power_EGID['EGID'].astype('int64')
    merged_gdf = gdf_cad.merge(power_EGID, left_on='RegBL_EGID', right_on='EGID', how='left')
    
    if not merged_gdf.empty:
        columns_to_add = ['Power[W]']
        for column in columns_to_add:
            gdf_cad[column] = merged_gdf[column]
    
    # Get the coordinates of the centrale
    central_node = (gdf_centrale['geometry'].x, gdf_centrale['geometry'].y)
    central_node_point = Point(*central_node)
    central_node_coordinates = [central_node_point.x, central_node_point.y]
    
    # Create an empty networkx graph
    graph = nx.Graph()
    
    # Add nodes to the graph
    for index, row in gdf_cad.iterrows():
        graph.add_node(index, geometry=row.geometry.centroid, power=row['Power[W]'], EGID=row['RegBL_EGID'], npid=index+1, Type='HX')
    
    # Add the imposed node to the graph
    imposed_node_index = len(gdf_cad)  # Index of the imposed node in the graph
    imposed_node_geometry = central_node_point  # Assumes a single imposed node
    graph.add_node(imposed_node_index, geometry=imposed_node_geometry, power=0, EGID=0, npid=imposed_node_index+1, Type='start heating station')
    
    # Get the coordinates of the geometries centroids
    coordinates = gdf_cad.geometry.centroid.apply(lambda point: (point.x, point.y)).tolist()
    
    # Calculate the pairwise distances between coordinates
    distances = distance_matrix(coordinates, coordinates)
    
    # Calculate the distances between the imposed node and all other geometries
    imposed_node_distances = [imposed_node_geometry.distance(geometry) for geometry in gdf_cad.geometry]
    distances = np.vstack([distances, imposed_node_distances])
    distances = np.hstack([distances, np.append(imposed_node_distances, 0).reshape(-1, 1)])
    
    # Create edges based on the minimum spanning tree algorithm
    minimum_spanning_tree = nx.minimum_spanning_tree(nx.from_numpy_array(distances))
    
    # Find the index of the closest neighbor to the imposed node
    closest_neighbor_index = None
    closest_neighbor_distance = float('inf')
    for node1, node2 in minimum_spanning_tree.edges:
        if node1 == imposed_node_index:
            distance = imposed_node_distances[node2]
            if distance < closest_neighbor_distance:
                closest_neighbor_index = node2
                closest_neighbor_distance = distance
        elif node2 == imposed_node_index:
            distance = imposed_node_distances[node1]
            if distance < closest_neighbor_distance:
                closest_neighbor_index = node1
                closest_neighbor_distance = distance
    
    # Perform a breadth-first search from the closest neighbor to assign key numbers to each node
    bfs_tree = nx.bfs_tree(minimum_spanning_tree, source=imposed_node_index)
    key_numbers = {node: 0 for node in bfs_tree.nodes}  # Initialize all key numbers to 0 for now
    parent_numbers = {node: 0 for node in bfs_tree.nodes}  # Initialize all parent numbers to 0 for now
    for node in bfs_tree.nodes:
        incoming_edges = list(bfs_tree.in_edges(node))  # Convert InEdgeDataView to a list
        if len(incoming_edges) > 0:
            parent_node = incoming_edges[0][0]  # Assumes a single incoming edge
            key_numbers[node] = key_numbers[parent_node] + 1
            parent_numbers[node] = parent_node
    
    # Assign the key numbers to the graph nodes
    nx.set_node_attributes(graph, key_numbers, name='key')
    nx.set_node_attributes(graph, parent_numbers, name='parent')
    
    # Sort the nodes based on their key numbers in descending order
    sorted_nodes = sorted(graph.nodes(data=True), key=lambda x: x[1]['key'], reverse=True)
    
    # Initialize a dictionary to store the cumulative power for each node
    cumulative_power = {}
    node_geometries = []
    node_coordinates_x = []
    node_coordinates_y = []
    node_npid = []
    node_EGID = []
    node_type = []
    
    # Iterate over the nodes in sorted order and calculate the cumulative power
    for node, data in sorted_nodes:
        cumulative_power[node] = data['power']
        node_npid.append(data['npid'])
        node_EGID.append(int(data['EGID']))
        node_type.append(data['Type'])
        node_geometries.append(data['geometry'])
        node_coordinates_x.append(round(data['geometry'].x,2))
        node_coordinates_y.append(round(data['geometry'].y,2))
         
    for node, data in sorted_nodes:
        node_parent = data['parent']
        if data['key'] != 0:
            cumulative_power[node_parent] += cumulative_power[node]
            
    # Create a new GeoDataFrame from the nodes geometries
    nodes_gdf = gpd.GeoDataFrame(geometry=node_geometries)
    nodes_gdf['coordinates_x'] = node_coordinates_x
    nodes_gdf['coordinates_y'] = node_coordinates_y
    nodes_gdf['npid'] = node_npid
    nodes_gdf['EGID'] = node_EGID
    nodes_gdf['Type'] = node_type
    
    # Iterate over the edges of the minimum spanning tree and size the lines based on cumulative power
    line_geometries = []
    line_powers = []
    line_nodes_start = []
    line_nodes_end = []
    line_pid = []
    line_length = []
    
    index = 0
    for edge in minimum_spanning_tree.edges:
        node1, node2 = edge
        key1 = graph.nodes[node1]['key']
        key2 = graph.nodes[node2]['key']
        npid1 = graph.nodes[node1]['npid']
        npid2 = graph.nodes[node2]['npid']
        geometry1 = graph.nodes[node1]['geometry']
        geometry2 = graph.nodes[node2]['geometry']
        line = LineString([geometry1, geometry2])
        line_power = cumulative_power[node1] if key1 > key2 else cumulative_power[node2]
        line_geometries.append(line)
        line_powers.append(line_power)
        line_nodes_start.append(npid1)
        line_nodes_end.append(npid2)
        line_pid.append(index)
        line_length.append(round(line.length,2))
        
        # Sizing of pipes
        for i in df_sizing.index:
            row = df_sizing.loc[i]
            power_kW = row['P [kW]']
            if line_power >= power_kW*1000:
                DN_index = i
        row = df_sizing.loc[DN_index]
        DN = row['DN']
        insul_th = row['t [mm]']/1000
        insul_k = row['U [W/mK]']
        
        index += 1
    
    # Create a new GeoDataFrame from the lines geometries
    lines_gdf = gpd.GeoDataFrame(geometry=line_geometries)
    lines_gdf['power_line[W]'] = line_powers
    lines_gdf['startpoint'] = line_nodes_start
    lines_gdf['endpoint'] = line_nodes_end
    lines_gdf['pid'] = line_pid
    lines_gdf['length[m]'] = line_length
    lines_gdf['DN'] = DN
    lines_gdf['insulation_thickness'] = insul_th
    lines_gdf['insulation_k_value'] = insul_k
    
    # Create a new DataFrame for the substations properties
    substations_df = gdf_cad[['RegBL_EGID','Power[W]']]
    
    # Create points and pipes for xml generation
    points = nodes_gdf.copy()
    pipes = lines_gdf.copy()
    substations = substations_df.copy()
    substations = substations.rename(columns={'RegBL_EGID':'EGID'})
    
    return graph, lines_gdf, nodes_gdf, points, pipes, substations

###################################################
# 
#              Plot Trace
#
##################################################

# graph, lines_gdf = get_trace(gdf_cad, gdf_centrale, power_EGID)

def plot_network_sizing(directory_path, graph, gdf_cad, lines_gdf, gdf_centrale):
    lines = lines_gdf.copy()
    lines['power_line[W]'] /= 1000
    lines = lines.rename(columns={'power_line[W]':'power_line[kW]'})

    fig, ax = plt.subplots(figsize=(10, 10), num='Network sizing')
    ax.set_axis_off()
    gdf_cad.plot(ax=ax, color='white', edgecolor='black', linewidth=1)
    ax.set_title("Network power sizing")
    lines.plot(ax=ax, column='power_line[kW]', cmap='autumn', linewidth=2, legend=True, legend_kwds={'label': "Power transferred through pipes [kW]"})
    gdf_centrale.plot(ax=ax, color='black', markersize=40)
    for node in graph.nodes:
        x, y = graph.nodes[node]['geometry'].coords[0]
        power_number = int(graph.nodes[node].get('power', '')/1000)
        text = ax.annotate(f"{power_number}kW", xy=(x, y), xytext=(0, 0), textcoords="offset points", color='black', fontsize=8)
        text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'))
    
    network_sizing_path = os.path.join(directory_path, 'Network_sizing.png')     
    plt.savefig(network_sizing_path, dpi=300, bbox_inches='tight')
    plt.show()

###################################################
# 
#              Plot Massflow
#
##################################################

def get_mass_flow(lines_gdf, scenario_path, sc_id):
    lines_mass_flow = lines_gdf.copy()
    pipes_results_path = os.path.join(scenario_path, f'Pipes_mass_flow_sc_{sc_id}.csv') 
    Pipes_mass_flow = pd.read_csv(pipes_results_path, delimiter=",")

    # Perform the merge on matching values "pid"
    merged_lines = lines_mass_flow.merge(Pipes_mass_flow, left_on='pid', right_on='pid', how='left')
    if not merged_lines.empty:
        column = ['Mass_flow[kg/s]']
        lines_mass_flow[column] = merged_lines[column]
    
    return lines_mass_flow

def get_supply_temp(nodes_gdf, scenario_path, sc_id):
    nodes_supply_temp = nodes_gdf.copy()
    points_results_path = os.path.join(scenario_path, f'Nodes_supply_temp_sc_{sc_id}.csv') 
    Points_supply_temp = pd.read_csv(points_results_path, delimiter=",")

    # Perform the merge on matching values "npid"
    merged_nodes = nodes_supply_temp.merge(Points_supply_temp, left_on='npid', right_on='npid', how='left')
    if not merged_nodes.empty:
        column = ['T_supply[°C]']
        nodes_supply_temp[column] = merged_nodes[column]

    return nodes_supply_temp

def Plot_Network_Data_v0(gdf_cad, lines_mass_flow, gdf_centrale):
    fig, ax = plt.subplots(figsize=(10, 10))
    gdf_cad.plot(ax=ax, color='white', edgecolor='black', linewidth=1)
    ax.set_axis_off()
    ax.set_title("Network pipes mass flow rate")
    lines_mass_flow.plot(ax=ax, column='Mass_flow[kg/s]', cmap='plasma', linewidth=2, legend=True, legend_kwds={'label': "Max Mass Flow rate [kg/s]"})
    gdf_centrale.plot(ax=ax, color='green', markersize=50)
    plt.show()

def plot_network_data(scenario_path, sc_id, graph, gdf_cad, 
                      lines_gdf, nodes_gdf, gdf_centrale, index_max):
    
    lines_mass_flow = get_mass_flow(lines_gdf, scenario_path, sc_id)
    nodes_supply_temp = get_supply_temp(nodes_gdf, scenario_path, sc_id)
    
    total_hours = int(index_max)  # Replace with the desired number of hours
    base_date = datetime.datetime(2023, 1, 1) #or any year that is not a leap year
    result_date_time = base_date + datetime.timedelta(hours=total_hours)
    result_date_time = result_date_time.strftime('%d.%m, %H:%M')
    
    fig, ax = plt.subplots(figsize=(10, 10), num=f'Network results scenario {sc_id}')
    gdf_cad.plot(ax=ax, color='white', edgecolor='black', linewidth=1)
    ax.set_axis_off()
    
    title = f"Network simulation results : scenario {sc_id}"
    # subtitle = f"Data obtained on {result_date_time}"
    print(f'scenario {sc_id} : {result_date_time}')
    ax.set_title(title, fontsize=12)
    # ax.text(0.5, 0.2, subtitle, transform=ax.transAxes, fontsize=10, ha='center', va='center')

    lines_mass_flow.plot(ax=ax, column='Mass_flow[kg/s]', cmap='winter', linewidth=2, legend=True, legend_kwds={'label': "Max mass flow rate [kg/s]"})
    gdf_centrale.plot(ax=ax, color='black', markersize=50)
 
    for node in graph.nodes:
        x, y = graph.nodes[node]['geometry'].coords[0]
        id_node = int(graph.nodes[node].get('npid', ''))
        T_supply = nodes_supply_temp[nodes_supply_temp['npid'] == id_node]['T_supply[°C]'].iloc[0]      
        text = ax.annotate(f"{T_supply:.1f}°C", xy=(x, y), xytext=(0, 0), textcoords="offset points", color='black', fontsize=8)
        text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'))

    network_results_path = os.path.join(scenario_path, 'Network_results.png')     
    plt.savefig(network_results_path, dpi=300, bbox_inches='tight')

    plt.show()
    
    
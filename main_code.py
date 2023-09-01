# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:16:28 2023

@author: Olivier Chavanne
"""

import geopandas as gpd
import pandas as pd
from shapely import box
import os
import matplotlib.pyplot as plt

# Local libraries
from enerCAD.building import generate_envelope
from enerCAD.building import generate_buildings
import enerCAD.xml as xml
import enerCAD.result as result
import enerCAD.network as network
import enerCAD.production as prod
import enerCAD.KPI as KPI
import enerCAD.IDC as IDC

# URL for RegBL API request
GEOADMIN_BASE_URL = "https://api.geo.admin.ch/rest/services/ech/MapServer/ch.bfs.gebaeude_wohnungs_register/"
    
##################################################
# 
#                  Functions
#
##################################################

def create_xml_root(xml_file_to_copy, climate_file, horizon_file):
    '''
    Parameters                                                          
    ----------
    xml_file_to_copy : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.

    Returns
    -------
    root : TYPE
        DESCRIPTION.
    district : TYPE
        DESCRIPTION.
    '''
    
    # Write XML file for CitySim :
    print("Writing XML file...")    
    # Add Root 
    root = xml.add_root()
    # Add Simulation days
    xml.add_simulation_days(root)
    # Add Climate
    xml.add_climate(root, climate_file)
    # Add District
    district = xml.add_district(root)
    
    # Horizon
    # read in the tab-separated file as a dataframe
    horizon_df = pd.read_csv(horizon_file, sep='\t', header=None)
    # assign column names to the dataframe
    horizon_df.columns = ['phi', 'theta']
    # Add Far field obstructions
    xml.add_far_field_obstructions(district, horizon_df)
    
    # Add all the composites and profiles, taken from a source XML
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Composite')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyYearProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DeviceType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'ActivityType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWYearProfile')
    
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Building')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DistrictEnergyCenter')
    
    print("Xml source copied")
    
    return root, district 

def Module_1(gpkg_filepath, GEOADMIN_BASE_URL,
             directory_path, xml_name,
             xml_base_file, climate_file, horizon_file,
             create_geometry_3D=False, calculate_volume_3D=False,
             EGID_column='RegBL_EGID'):
    '''
    Parameters
    ----------
    gpkg_filepath : TYPE
        DESCRIPTION.
    GEOADMIN_BASE_URL : TYPE
        DESCRIPTION.
    directory_path : TYPE
        DESCRIPTION.
    xml_file_to_create : TYPE
        DESCRIPTION.
    xml_base_file : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.
    create_geometry_3D : TYPE, optional
        DESCRIPTION. The default is False.
    calculate_volume_3D : TYPE, optional
        DESCRIPTION. The default is False.
    EGID_column : TYPE, optional
        DESCRIPTION. The default is 'RegBL_EGID'.

    Returns
    -------
    None.
    '''
    
    ### Exctract geopackage ###
    
    print("Exctracting geopackage layers...")
    
    # MO Cadaster
    MO_all = gpd.read_file(gpkg_filepath, layer = "zone_tout")
    MO_dhn = gpd.read_file(gpkg_filepath, layer = "zone_cad")
    centrale = gpd.read_file(gpkg_filepath, layer = "centrale")
    EGID_column = 'RegBL_EGID'
    
    # Split Multipolygons into Polygons
    zone_all = MO_all.explode(index_parts=False)
    zone_dhn = MO_dhn.explode(index_parts=False)
    
    # List containing EGID of buildings to simulate
    EGID_list = MO_dhn[EGID_column].tolist()
    
    # Save EGID list of buildings connected to CAD
    df_EGID = pd.DataFrame(EGID_list)
    df_EGID.columns = ['EGID']
    EGID_path = os.path.join(directory_path, 'EGID.csv')     
    df_EGID.to_csv(EGID_path, index=False)
    print("EGID.csv created")
    
    # Swissbuildings3D
    print("Swissbuildings3D processing...")
    try:
        floor_data = gpd.read_file(gpkg_filepath, layer = "floor")
        roof_data = gpd.read_file(gpkg_filepath, layer = "roof")
        wall_data = gpd.read_file(gpkg_filepath, layer = "wall")
        
        # Filter on the zone with 10m buffer around surrounding square box 
        zone_bounds = MO_all.geometry.buffer(10).values.total_bounds
        zone_box = box(zone_bounds[0], zone_bounds[1], zone_bounds[2], zone_bounds[3])
        
        # Cut swissbuildings3D to zone of concern
        floor_data_intersection = floor_data[floor_data.geometry.intersects(zone_box)]
        roof_data_intersection = roof_data[roof_data.geometry.intersects(zone_box)]
        wall_data_intersection = wall_data[wall_data.geometry.intersects(zone_box)]
    
        # Split Multipolygons into Polygons
        zone_floor = floor_data_intersection.explode(index_parts=True).reset_index()
        zone_roof = roof_data_intersection.explode(index_parts=True).reset_index()
        zone_wall = wall_data_intersection.explode(index_parts=True).reset_index()
        print('Swissbuildings3D cut to zone of interest \n')
    
    except: print('Error : Swissbuildings3D not provided')

    ### Envelope processing ###
    
    try:
        # Get z coordinates of 1st vertex from 1st surface of 1st building's floor polygon as altitude by default for MO footprints
        altitude_default = zone_floor.loc[0].geometry.exterior.coords[0][2]
    except:
        altitude_default = 0
    
    # Create DataFrames containing all necessary information for each building
    print("Creating Buildings GeoDataFrame...")
    footprints, buildings = generate_buildings(zone_all, EGID_list, GEOADMIN_BASE_URL, altitude_default,
                                               create_geometry_3D, calculate_volume_3D, zone_floor, zone_roof, zone_wall)
    print("Buildings GeoDataFrame created \n") 
    
    # Generate the envelope surfaces
    print("Generating Buildings envelope...")
    envelope, buildings_volume_3D, center_coordinates = generate_envelope(footprints, buildings, calculate_volume_3D)
    print("Envelope created \n")
    
    # Merge "volume_3D" and "n_occupants" to main buildings geodataframe according to 'bid'
    merged_buildings = buildings.merge(buildings_volume_3D, left_on='bid', right_on='bid', how='left')    
    if not merged_buildings.empty:
        columns_to_add = ['volume_3D', 'n_occupants']
        for column in columns_to_add:
            buildings[column] = merged_buildings[column]
        print("Buildings 3D volume calculated and merged \n")
    
    ### Buildings XML processing ###
        
    root, district = create_xml_root(xml_base_file, climate_file, horizon_file)
    
    print("Adding buildings...")
    # Add the buildings
    xml.add_all_buildings(district, buildings, envelope, center_coordinates)
        
    # Write XML file
    xml_path = os.path.join(directory_path, xml_name+".xml")     
    xml.write_xml_file(root, xml_path)
    print(f"{xml_name}.xml file created \n")
    
    return buildings, zone_dhn, centrale
   
def simulate_citysim(directory_path, xml_file, citysim_filepath):
    '''
    Parameters
    ----------
    xml_file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    '''
    
    import subprocess
    import time
    start = time.time()
    print('Process started')
    print(f'Simulation of {xml_file}.xml...')

    #run CitySim.exe with xml file
    xml_path = os.path.join(directory_path, xml_file+".xml")
    result = subprocess.run([citysim_filepath, '-q', f"{xml_path}"])
    
    end = time.time()
    duration = end - start
    m, s = divmod(duration, 60)
    print('Simulation ended. Time :', "%.0f" %m,'min', "%.0f" %s,'s \n')
    

#------------------Part 2 iterating----------------------------------------------------------

def Module_2(directory_path, xml_name, zone_dhn, centrale,
             xml_DHN, climate_file, horizon_file,
             scenarios_list):   

    # Open the 1st iteration results file (without DHN)
    results_filepath = os.path.join(directory_path, xml_name+"_TH.out")
    results = pd.read_csv(results_filepath, delimiter="\t")
    results = results.set_index("#timeStep")
    print(f'{xml_name}_TH.out opened')
    
    # Get P_max for every building in the network
    power_EGID = result.get_Pmax_per_EGID(results)
    power_EGID_path = os.path.join(directory_path, 'Power_EGID.csv')
    power_EGID.to_csv(power_EGID_path, index=False)
        
    # Create and size network
    graph, lines_gdf, nodes_gdf, points, pipes, substations = network.get_trace(zone_dhn, centrale, power_EGID)
    
    # Save to csv file
    points_path = os.path.join(directory_path, 'Points.csv')
    pipes_path = os.path.join(directory_path, 'Pipes.csv')
    substations_path = os.path.join(directory_path, 'Substations.csv')
    points.to_csv(points_path, index=False)
    pipes.to_csv(pipes_path, index=False)
    substations.to_csv(substations_path, index=False)
    print("csv files saved")
    
    # Get Load duration curve of the network
    load_curve = result.get_thermal_load_curve(results)
    load_curve_path = os.path.join(directory_path, 'Load_curve.csv')
    load_curve.to_csv(load_curve_path, index=False)
     
    # Compute production scenarios
    scenarios = prod.get_scenarios(load_curve)
    
    # Calculate production water storage
    volume_storage, capacity_storage = prod.get_storage(load_curve)
    
    for i in range(len(scenarios)):
        scenario = scenarios[i]
        sc_id = scenario[0]
        if sc_id in scenarios_list:

            ### Scenarios XML processing ###
            
            xml_to_copy_path = os.path.join(directory_path, xml_name+'.xml' )
            root, district = create_xml_root(xml_to_copy_path, climate_file, horizon_file)
               
            # Add District heating center
            district_heating_center = xml.add_district_heating_center(district)
            
            # Add Network (nodes and pipes)
            xml.add_network(district_heating_center, points=points.copy(), pipes=pipes.copy())

            # Change Boilers into Substations
            xml.change_boiler_to_substation(district, substations=substations.copy(), points=points.copy())
            
            # Add Thermal station
            ts_node_id = points.loc[points['Type']=='start heating station']['npid'].iloc[0]
            network_p_max = pipes['power_line[W]'].max()
            production_stages = [scenario[1],scenario[2]]
            
            # Efficiency data
            technology_parameters = pd.read_csv("KPI.csv", delimiter=",")
            eff_columns = ['T_eff','efficiency']
            efficiency_parameters = technology_parameters[eff_columns]
            
            xml.add_thermal_station(district_heating_center, ts_node_id, p_max=network_p_max, 
                                    c_storage=capacity_storage, efficiencies=efficiency_parameters, stages=production_stages)

            # Write XML file
            scenario_path = os.path.join(directory_path,f"Scenario_{sc_id}")
            os.makedirs(scenario_path, exist_ok=True)
            xml_to_create_path = os.path.join(scenario_path, xml_DHN+f"_sc_{sc_id}"+".xml")
            xml.write_xml_file(root, xml_to_create_path)
            print(f'{xml_DHN}_sc_{sc_id}.xml file created \n')

    return graph, lines_gdf, nodes_gdf, results, scenarios, volume_storage

#--------------------- Part 3

def Module_results_network(scenario_path, sc_id, xml_DHN, zone_dhn, centrale, graph, lines_gdf, nodes_gdf):
    
    # Open the 2nd iteration results file (with DHN)
    results_filepath = os.path.join(scenario_path, xml_DHN+f"_sc_{sc_id}"+"_TH.out")
    results = pd.read_csv(results_filepath, delimiter="\t")
    results_final = results.set_index("#timeStep")
        
    # Get Data for every node and pipe in the network
    Pipes_mass_flow, Nodes_supply_temp, index_max = result.get_network_data_max(results_final)
    
    # Save to csv file
    Pipes_mass_flow_filepath = os.path.join(scenario_path, f'Pipes_mass_flow_sc_{sc_id}.csv')
    Nodes_supply_temp_filepath = os.path.join(scenario_path, f'Nodes_supply_temp_sc_{sc_id}.csv')
    Pipes_mass_flow.to_csv(Pipes_mass_flow_filepath, index=False) 
    Nodes_supply_temp.to_csv(Nodes_supply_temp_filepath, index=False)

    return results_final, Pipes_mass_flow, Nodes_supply_temp, index_max

#--------------------- KPI calculation

def Module_KPI(results_production, volume_storage, 
               scenarios, sc_id, scenario_path, do_plot):

    # Get production consumption in kWh
    pump_cons, fuel_cons, elec_cons, th_prod, df_th_prod, df_elec = result.get_energy_data(results_production, sc_id, scenarios)

    # Save thermal results to csv file
    th_production_results_filepath = os.path.join(scenario_path, f'Thermal_sc_{sc_id}.csv')
    df_th_prod.to_csv(th_production_results_filepath, index=True) 

    if do_plot == True:
    # Plot energy production data
        print('Energy production plot...')
        result.plot_energy_data(results_production, sc_id, scenarios, scenario_path)

    # Calculate KPI (key performance indicators)
    technology_parameters = pd.read_csv("KPI.csv", delimiter=",")

    df_KPI = KPI.calculate_KPI(sc_id, scenarios, volume_storage, technology_parameters,
                                                 pump_cons, fuel_cons, elec_cons, th_prod)
    print('KPI calculated')
    
    # Save electrical results to csv file
    electricity_results_filepath = os.path.join(scenario_path, f'Electricity_sc_{sc_id}.csv')
    df_elec.to_csv(electricity_results_filepath, index=True) 

    # Save KPI results to csv file
    KPI_results_filepath = os.path.join(scenario_path, f'KPI_results_sc_{sc_id}.csv')
    df_KPI.to_csv(KPI_results_filepath, index=False) 
    
    return df_KPI


##################################################
# 
#         Information to provide
#
##################################################

# Geopackage filepath
gpkg_filepath = r"Orbe.gpkg"                                   #TODO
# gpkg_filepath = r"zone_Delémont.gpkg"                                   #TODO
# gpkg_filepath = r"test_multi_prod.gpkg"                                   #TODO

# Create geometry with swissbuildings3D
create_geometry_3D = False                                               #TODO

# Calculate volume from swissbuildings3D
calculate_volume_3D = False                                               #TODO

# CitySim.exe filepath
citysim_filepath = "C:/Users/Olivier Chavanne/Documents/CitySim_solver/v7/bin/CitySim.exe" #TODO

# XML name to export
directory_path = "Orbe_2.5D"

os.makedirs(directory_path, exist_ok=True)
                                      
xml_name = directory_path                                           #TODO
xml_DHN = "DHN_"+xml_name

# XML source files
xml_base_file = 'xml_base.xml'                                          #TODO
climate_file = 'Bex_2030.cli'                                           #TODO
horizon_file = 'Bex.hor'                                                #TODO

# Scenarios to simulate
scenarios_list = [1,2,3,4,5,6]                              #TODO

# For IDC calculation
measures_filepath = r"Regener_2023-05-03_Orbe_Pully_Bex.xlsx" #TODO

do_plot = True

def main(): 
    
    # Generate individual buildings XML
    print('***Module 1*** \n')
    buildings, zone_dhn, centrale = Module_1(gpkg_filepath, GEOADMIN_BASE_URL, 
                                             directory_path, xml_name,
                                             xml_base_file, climate_file, horizon_file,
                                             create_geometry_3D, calculate_volume_3D,
                                             EGID_column='RegBL_EGID')
 
    # 1st CitySim simulation
    simulate_citysim(directory_path, xml_name, citysim_filepath)
    
    # Generate DHN XML for each scenario
    print('***Module 2*** \n')
    graph, lines_gdf, nodes_gdf, results, scenarios, volume_storage = Module_2(directory_path, xml_name, 
                                                                        zone_dhn, centrale, 
                                                                        xml_DHN, climate_file, horizon_file,
                                                                        scenarios_list)
    # Intermediate results
    network.plot_network_sizing(directory_path, graph, zone_dhn, lines_gdf, centrale)     
    result.plot_load_curve(results, directory_path)
    
    KPI_result_list = []
        
    # DHN simulation for each scenario
    for i in range(len(scenarios_list)):
        sc_id = scenarios_list[i]
        print(f'***Scenario {sc_id}*** \n')
        
        # CitySim simulation        
        scenario_path = os.path.join(directory_path,f"Scenario_{sc_id}")
        simulate_citysim(scenario_path, f'{xml_DHN}_sc_{sc_id}', citysim_filepath)
                
        # Analysing final results  
        print('Results processing...')
        results_prod_path = os.path.join(scenario_path,  f'{xml_DHN}_sc_{sc_id}_TH.out')
        results_production = pd.read_csv(results_prod_path, delimiter="\t")
            
        # KPI calculation
        df_KPI = Module_KPI(results_production, volume_storage, 
                            scenarios, sc_id, scenario_path, do_plot)
        
        KPI_result_list.append(df_KPI)

        # Network final plots
        results_final, Pipes_mass_flow, Nodes_supply_temp, index_max = Module_results_network(scenario_path, sc_id, xml_DHN, 
                                                                                zone_dhn, centrale, 
                                                                                graph, lines_gdf, nodes_gdf)   
        if do_plot == True:
            network.plot_network_data(scenario_path, sc_id, graph, zone_dhn, 
                                      lines_gdf, nodes_gdf, centrale, index_max)
        
        plt.close("all")
        print(f"Scenario {sc_id} processed \n")
    
    if len(KPI_result_list)>1:
        KPI.plot_KPI(KPI_result_list, scenarios_list, directory_path)
    
    print("***Overall processing finished***")
    print(f"Find all results and graphs in directory : {directory_path}")

if __name__ == "__main__":
    plt.close("all")
    
    import subprocess
    import time
    start_overall = time.time()
    print('Main code started')
    print('-----------------')
    
    main()    

    print('-----------------')
    print('Main code ended')
    print('-----------------')
    end_overall = time.time()
    duration_overall = end_overall - start_overall
    m, s = divmod(duration_overall, 60)
    print('Overall run time :', "%.0f" %m,'min', "%.0f" %s,'s \n')















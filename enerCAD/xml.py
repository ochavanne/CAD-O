# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 13:18:03 2021

@author: Roberto
"""
from xml.etree.ElementTree import Element, SubElement, tostring, parse
from xml.dom.minidom import parseString
import numpy as np

'''
TO DO:
    
    - thank Pierre
    - unify egid, ssid

'''

def prettify(element):
    """
    Return a pretty-printed XML string for the Element.
    """
    rough_string = tostring(element, 'utf-8')
    reparsed = parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def write_xml_file(root, filename):
    text = prettify(root)
    with open(filename, "w") as f:
        f.write(text)
    # Remove empty lines.
    filedata = ""
    with open(filename, "r") as infile:
        for line in infile.readlines():
            if line.strip():  # if striped line is non empty
                filedata += line
    with open(filename, "w") as outfile:
        outfile.write(filedata)

def get_district_child_from_xml(xml_file, tag):
    tree = parse(xml_file)
    root = tree.getroot()
    file_district = root.find('District')
    elements = file_district.findall(tag)
    return elements


def add_child_from_xml_to_district(district, xml_file, tag):
    elements = get_district_child_from_xml(xml_file, tag)
    for element in elements:
        district.append(element)


def add_composites_from_xml(district, xml_file):
    add_child_from_xml_to_district(district, xml_file, 'Composite')


def add_far_fields_obstructions_from_xml(district, xml_file):
    add_child_from_xml_to_district(district, xml_file, 'FarFieldObstructions')


def add_profiles_from_xml(district, xml_file):
    add_child_from_xml_to_district(district, xml_file, 'Composite')
    add_child_from_xml_to_district(district, xml_file, 'OccupancyDayProfile')
    add_child_from_xml_to_district(district, xml_file, 'OccupancyYearProfile')
    add_child_from_xml_to_district(district, xml_file, 'DeviceType')
    add_child_from_xml_to_district(district, xml_file, 'ActivityType')
    add_child_from_xml_to_district(district, xml_file, 'DHWDayProfile')
    add_child_from_xml_to_district(district, xml_file, 'DHWYearProfile')


def add_root():
    citysim = Element("CitySim")
    return citysim


def add_simulation_days(citysim, begin_month=1, begin_day=1,
                        end_month=12, end_day=31):

    dict_simulation = {"beginMonth": str(begin_month),
                       "beginDay": str(begin_day),
                       "endMonth": str(end_month),
                       "endDay": str(end_day)}
    simulation = SubElement(citysim, "Simulation", dict_simulation)
    return simulation


def add_climate(citysim, filename):
    dict_climate = {"location": filename}
    climate = SubElement(citysim, "Climate", dict_climate)


def add_district(citysim):
    district = SubElement(citysim, "District")
    return district

def add_far_field_obstructions(district, horizon_df):
    far_field_obstructions = SubElement(district, "FarFieldObstructions")
    for i in horizon_df.index:
        row = horizon_df.loc[i]
        phi = row['phi']
        theta = row['theta']
        dict_horizon = {"phi": str(phi), "theta": str(theta)}
        horizon = SubElement(far_field_obstructions, "Point", dict_horizon)


def add_composite(district, composite_id, composite_name, layers_df):
    dict_composite = {"id": str(composite_id), "name": composite_name}
    composite = SubElement(district, "Composite", dict_composite)
    for i in layers_df.index:
        row = layers_df.loc[i]
        name = row['name']
        thickness = row['thickness']
        conductivity = row['conductivity']
        cp = row['cp']
        density = row['density']
        dict_layer = {"Name": name, "Thickness": str(thickness),
                      "Conductivity": str(conductivity), "Cp": str(cp),
                      "Density": str(density)}
        layer = SubElement(composite, "Layer", dict_layer)


def add_building(district, row, volume, simulate=True, tmin=21, tmax=26, 
                 blinds_lambda=0.2,blinds_irradiance_cutoff=150):
    # Building characteristics
    egid = row['egid']
    bid = row['bid']
    ventilation_rate = row['Ninf']

    dict_building = {"id": str(bid), "key": str(egid), "Ninf": str(ventilation_rate),
                     "Vi": str(volume), "Tmin": str(tmin), "Tmax": str(tmax),
                     "BlindsLambda": str(blinds_lambda),
                     "BlindsIrradianceCutOff": str(blinds_irradiance_cutoff),
                     "Simulate": str(simulate).lower()}

    building = SubElement(district, "Building", dict_building)
    return building


def add_heat_tank(building, v=50*1e-3, phi=200, rho=1000, cp=4180, tmin=20,
                  tmax=35, tcrit=90):

    dict_tank = {"V": str(v), "phi": str(phi), "rho": str(rho),
                 "Cp": str(cp), "Tmin": str(tmin), "Tmax": str(tmax),
                 "Tcritical": str(tcrit)}

    heat_tank = SubElement(building, "HeatTank", dict_tank)
    return heat_tank


def add_dhw_tank(building, v=0.2, phi=2.5, rho=1000, cp=4180, tmin=60,
                 tmax=70, tcrit=90, tinlet=5):

    dict_tank = {"V": str(v), "phi": str(phi), "rho": str(rho),
                 "Cp": str(cp), "Tmin": str(tmin), "Tmax": str(tmax),
                 "Tcritical": str(tcrit), "Tinlet": str(tinlet)}

    dhw_tank = SubElement(building, "DHWTank", dict_tank)
    return dhw_tank


def add_cool_tank(building, v=20, phi=20, rho=1000, cp=4180, tmin=5,
                  tmax=20):

    dict_tank = {"V": str(v), "phi": str(phi), "rho": str(rho),
                 "Cp": str(cp), "Tmin": str(tmin), "Tmax": str(tmax)}

    cool_tank = SubElement(building, "CoolTank", dict_tank)
    return cool_tank


def add_heat_source(building, begin_day=1, end_day=365):
    '''
    begin_day : The default is January 1st.
    end_day : The default is December 31st.
    '''
    
    dict_heat_source = {"beginDay": str(begin_day), "endDay": str(end_day)}
    heat_source = SubElement(building, "HeatSource", dict_heat_source)
    return heat_source


def add_boiler(heat_source, pmax=50000, eta_th=0.96):
    dict_boiler = {"name": "Boiler", "Pmax": str(pmax), "eta_th": str(eta_th)}
    boiler = SubElement(heat_source, "Boiler", dict_boiler)
    return boiler


def add_zone(building, net_volume, zone_id=0, psi=0.2, ground_floor=True):

    dict_zone = {"id": str(zone_id), "volume": str(net_volume),
                 "Psi": str(psi), "groundFloor": str(ground_floor).lower()}

    zone = SubElement(building, "Zone", dict_zone)
    return zone


def add_imposed_heat_demand(building, values, start_day=1, start_hour=1,
                            end_day=365, end_hour=24):

    count = 0
    demand_dict = {}
    for i in range(start_day, end_day + 1):
        if start_day == end_day:
            daily_hours_range = range(start_hour, end_hour + 1)
        elif i == start_day and start_day != end_day:
            daily_hours_range = range(start_hour, 25)
        elif i == end_day and start_day != end_day:
            daily_hours_range = range(1, end_hour + 1)
        else:
            daily_hours_range = range(1, 25)
        for j in daily_hours_range:
            key = 'd{}h{}'.format(i, j)
            value = values.loc[count]
            demand_dict[key] = str(value)
            count += 1
    heat_demand = SubElement(building, "ImposedHeatDemand", demand_dict)
    return heat_demand


def add_occupants(zone, number_of_occupants=5, building_type=1, activity_type=None,
                  dhw_type=None, stochastic=False):

    dict_occupants = {"n": str(number_of_occupants),
                      "type": str(building_type), "Stochastic": str(stochastic).lower(),
                      # "activityType":  str(activity_type),
                      # "DHWType": str(dhw_type)
                      }
    if activity_type:
        dict_occupants["activityType"] = str(activity_type)
    if dhw_type:
        dict_occupants["DHWType"] = str(dhw_type)
    occupants = SubElement(zone, "Occupants", dict_occupants)
    return occupants

def add_surfaces(zone, envelope, center_coordinates):
    x_center = center_coordinates[0]
    y_center = center_coordinates[1]
    for r in envelope.index:
        row = envelope.loc[r]
        geometry = row['geometry']
        class_id = row['class_id']
        glazing_ratio = row["glazing_ratio"]
        glazing_u_value = row["glazing_u_value"]
        glazing_g_value = row["glazing_g_value"]
        openable_ratio = row["openable_ratio"]
        shortwave_reflectance = row["shortwave_reflectance"]
        surface_type = int(row["surface_type"])

        if class_id == 33:
            dict_surface = {"id": str(r), "type": str(surface_type)}
            surface = SubElement(zone, "Floor", dict_surface)

        elif class_id == 34:
            dict_surface = {"id": str(r), "type": str(surface_type),
                            "ShortWaveReflectance": str(shortwave_reflectance),
                            "GlazingRatio": str(glazing_ratio),
                            "GlazingUValue": str(glazing_u_value),
                            "GlazingGValue": str(glazing_g_value),
                            "OpenableRatio": str(openable_ratio)}
            surface = SubElement(zone, "Wall", dict_surface)

        elif class_id == 35:
            dict_surface = {"id": str(r), "type": str(surface_type  ),
                            "ShortWaveReflectance": str(shortwave_reflectance),
                            "GlazingRatio": str(glazing_ratio),
                            "GlazingUValue": str(glazing_u_value),
                            "GlazingGValue": str(glazing_g_value),
                            "OpenableRatio": str(openable_ratio)}
            surface = SubElement(zone, "Roof", dict_surface)
            
        else:
            raise ValueError("Surface class not understood.")

        # Add points translated to center of the scene
        for p in range(len(geometry.exterior.coords)-1):
            point_name = "V{}".format(p)
            coordinates = geometry.exterior.coords[p]
            x = str(coordinates[0]-x_center)
            y = str(coordinates[1]-y_center)
            z = str(coordinates[2])

            dict_point = {"x": x, "y": y, "z": z}
            point = SubElement(surface, point_name, dict_point)

def add_ground(district, terrain_df, groundtype=37, detailedSimulation=False, ShortWaveReflectance=0.3):
    groundsurface = SubElement(district, "GroundSurface")
    for r in terrain_df.index:
        geometry = terrain_df["geometry"].loc[r]
        dict_surface = {"id": str(r),
                        "ShortWaveReflectance":str(ShortWaveReflectance),
                        "type":str(groundtype),
                        "detailedSimulation":str(detailedSimulation).lower()}

        surface = SubElement(groundsurface, "Ground", dict_surface)
        # Add points
        for p in range(len(geometry.exterior.coords)-1):
            point_name = "V{}".format(p)
            coordinates = geometry.exterior.coords[p]
            x = str(coordinates[0])
            y = str(coordinates[1])
            z = str(coordinates[2])

            dict_point = {"x": x, "y": y, "z": z}
            point = SubElement(surface, point_name, dict_point)

def add_all_buildings(district, buildings, envelope, center_coordinates=(0,0)):
    for i in buildings.index:
        row = buildings.loc[i]
        
        # Get volume from Swissbuildings3D if provided
        volume_MO = row['volume_MO']
        volume_3D = row['volume_3D']
        if volume_3D == 0:
            volume = volume_MO
        else:
            volume = volume_3D
            
        # Add building with according simulation status
        Simulate_status = row['Simulate_status']
        tmin = row['Tmin']
        if Simulate_status == True:
            building = add_building(district, row, volume, tmin=tmin, tmax=26, blinds_lambda=0.2,
                     blinds_irradiance_cutoff=150)
        else:
            building = add_building(district, row, volume, simulate=False)
        zone = add_zone(building, volume)
        
        # Activity profile according to building type (no profile for 11 sports installations)
        building_type = row['building_type']
        if building_type == 11:
            activity_type = None
        else:
            activity_type = building_type
        
        # DHW profile according to building type
        dhw_type = row['building_type']

        # Add heat tank with temperature setpoint of heat supplier depending on year of construction
        year = row['year']
        # Radiator
        if year < 1990:
            add_heat_tank(building, tmin=50, tmax=60, tcrit=90) #TODO
        # Underfloor heating
        else:
            add_heat_tank(building, tmin=35, tmax=40, tcrit=90)        
        
        # Add DHW tank according to number of occupants (0.05 m3/person, 3 m3 max.)
        n_occupants = row['n_occupants']
        if n_occupants != 0:
            vol_DHW = 50*1e-3*n_occupants
            if vol_DHW > 3:
                vol_DHW = 3   
            add_dhw_tank(building, v=vol_DHW)
        add_occupants(zone, n_occupants, building_type, activity_type, dhw_type, stochastic=False)
        
        # Add boiler as heat source of 10 MW
        heat_source = add_heat_source(building)
        add_boiler(heat_source, pmax=10e6)

        # Add building's envelope surfaces 
        bid = row['bid']
        e = envelope[envelope['bid'] == bid]
        add_surfaces(zone, e, center_coordinates)
        


# DHN ########################################################################

def add_district_heating_center(district, cp_water=4180, rho=990,
                                mu=0.0004):

    d = {"id": str(0), "Cp": str(cp_water), "rho": str(rho), "mu": str(mu)}
    district_heating_center = SubElement(district, 'DistrictEnergyCenter', d)
    return district_heating_center

def add_thermal_station(district_heating_center, node_id,
                        rho=999, n0=5000, a0_n0=1247180, a1_n0=-1640.236,
                        a2_n0=-0.00016031, e_pump=0.6,  begin_day=1, end_day=365,
                        low_supply_temp=90, high_supply_temp=75,
                        low_ext_temp_limit=-10, high_ext_temp_limit=15,
                        high_temperature_drop=20, start_summer_threshold=16,
                        end_summer_threshold=5, p_max=900000, delta_p=200000,
                        dp_type='constant', storage_type='seasonalStorageHeating',
                        kvMax=1000, temp_storage=75, c_storage=1e7,
                        efficiencies=None, stages=None):
    
    print("Preparing Thermal Station...")
    
    # Set thermal station parameters
    d_thermal_station = {"linkedNodeId": str(node_id), "beginDay": str(begin_day), 
                         "endDay": str(end_day), "type": str(storage_type), "kvMax": str(kvMax)}
    thermal_station = SubElement(district_heating_center, "ThermalStation",
                                 d_thermal_station)
    
    # Set temperature parameters
    d_temp_setpoint = {"type": "affineWinterConstantSummer",    
                        "lowExtTemp": str(low_ext_temp_limit),
                        "highExtTemp": str(high_ext_temp_limit),
                        "lowExtTempSupplyTemp": str(low_supply_temp),
                        "highExtTempSupplyTemp": str(high_supply_temp),
                        "startSummerTempThreshold": str(start_summer_threshold),
                        "endSummerTempThreshold": str(end_summer_threshold)}
    temp_setpoint = SubElement(thermal_station, "TemperatureSetpoint",
                                d_temp_setpoint)
    
    # Set pressure parameters
    if dp_type == 'affine': 
        mass_flows = np.array([5, 10, 64, 118])*rho/3600
        pressure_diffs = [180000, 250000, 390000, 480000]
        d_pressure_setpoint = {"type": "affine", "massFlow1": str(mass_flows[0]),
                                "pressureDiff1": str(pressure_diffs[0]),
                                "massFlow2": str(mass_flows[1]),
                                "pressureDiff2": str(pressure_diffs[1]),
                                "massFlow3": str(mass_flows[2]),
                                "pressureDiff3": str(pressure_diffs[2]),
                                "massFlow4": str(mass_flows[3]),
                                "pressureDiff4": str(pressure_diffs[3])}
        pressure_setpoint = SubElement(thermal_station, "PressureSetpoint",
                                        d_pressure_setpoint)
    elif dp_type == 'constant':
        d_pressure_setpoint = {"type": "constant",
                                "targetPressureDiff": str(delta_p)}
        pressure_setpoint = SubElement(thermal_station, "PressureSetpoint",
                                        d_pressure_setpoint)
    else:
        raise ValueError('dp_type must be either "constant" or "affine"')

    # Set pump parameters #TODO
    d_pump = {"n0": str(n0), "a0": str(a0_n0),
              "a1": str(a1_n0), "a2": str(a2_n0)}
    ts_pump = SubElement(thermal_station, "Pump", d_pump)

    d_epump = {"type": "constant", "efficiencyPump": str(e_pump)}
    epump = SubElement(ts_pump, "EfficiencyPump", d_epump)

    # Set storage parameters #TODO
    d_storage = {"type": "simple", "initialTemperature": str(high_supply_temp),
                  "heatCapacity": str(c_storage)}
    ts_storage = SubElement(thermal_station, "Storage", d_storage)

    # d_boiler = {"Pmax": str(p_max/2), "eta_th": "0.95"}
    # ts_boiler = SubElement(thermal_station, "Boiler", d_boiler)
    
    # d_chp = {"Pmax": str(p_max/2), "eta_th": "0.35", "eta_el": "0.6", "minPartLoadCoeff": "0.2"}
    # ts_chp = SubElement(thermal_station, "CHP", d_chp)

    # Add heat production             
    if stages == None:
        d_boiler = {"Pmax": str(p_max), "eta_th": "0.95"}
        ts_boiler = SubElement(thermal_station, "Boiler", d_boiler)
    
    else:
        # Set production units parameters
        for i in range(len(stages[0])):
            unit_power = stages[0][i]
            unit_type = stages[1][i]
            
            if unit_type == "CHP":
                eta_th = efficiencies[efficiencies['T_eff']=='CHP_th']['efficiency'].iloc[0]
                eta_el = efficiencies[efficiencies['T_eff']=='CHP_el']['efficiency'].iloc[0]
                d_CHPHP = {"Pmax": str(unit_power), "eta_th": str(eta_th), "eta_el": str(eta_el), "minPartLoadCoeff":"0.0"} #TODO
                ts_CHPHP = SubElement(thermal_station, "CHP", d_CHPHP)
    
            elif unit_type == "Heat_Pump_Water":
                eta_tech = efficiencies[efficiencies['T_eff']=='Heat_Pump_Water_eta_tech']['efficiency'].iloc[0]
                COP = efficiencies[efficiencies['T_eff']=='Heat_Pump_Water_COP']['efficiency'].iloc[0]
                P_el = unit_power/COP
                d_HeatPump = {"Pmax": str(P_el), "eta_tech": str(eta_tech), "Ttarget":"80", "Tsource":"20"} 
                ts_HeatPump = SubElement(thermal_station, "HeatPump", d_HeatPump)
        
            elif unit_type == "Wood_boiler":
                eta_th = efficiencies[efficiencies['T_eff']=='Wood_boiler']['efficiency'].iloc[0]
                d_WoodBoiler = {"Pmax": str(unit_power), "eta_th": str(eta_th), "name":"Wood Boiler"} #TODO
                ts_WoodBoiler = SubElement(thermal_station, "Boiler", d_WoodBoiler)   
    
            elif unit_type == "Gas_boiler":
                eta_th = efficiencies[efficiencies['T_eff']=='Gas_boiler']['efficiency'].iloc[0]
                d_GasBoiler = {"Pmax": str(unit_power), "eta_th": str(eta_th),"name":"Gas Boiler"} #TODO
                ts_GasBoiler = SubElement(thermal_station, "Boiler", d_GasBoiler)  


def add_network(district_heating_center, points, pipes):
    '''

    '''
    print("Preparing Network...")

    network = SubElement(district_heating_center,
                         "Network", {"soilkValue": "0.5"})

    # Set nodes parameters
    points = points.sort_values(['npid']).reset_index(drop=True)
    for index_node in points.index:
        current_row = points.loc[index_node]
        node_id = current_row['npid']
        node_coordinates_x = current_row['coordinates_x']
        node_coordinates_y = current_row['coordinates_y']
        node_EGID = current_row['EGID']
        node_type = current_row['Type']
        # Node connected to Thermal station
        if node_type == 'start heating station':
            pair_type = "ThermalStationNodePair"
        # Node connected to Substation
        elif node_type == 'HX':
            pair_type = "SubstationNodePair"
        # Node connecting pipes without heat exchange
        else:
            pair_type = "NodePair"
        d_point = {"id": str(node_id), "key": str(node_EGID), 
                   "x": str(node_coordinates_x), "y": str(node_coordinates_y), "z": str(0)}
        point_pair = SubElement(network, pair_type, d_point)

    # Set pipes parameters
    for index_pipe in pipes.index:
        current_row = pipes.loc[index_pipe]
        pipe_start_point = current_row['startpoint']
        pipe_end_point = current_row['endpoint']
        pipe_id = current_row['pid']
        pipe_length = current_row['length[m]']
        pipe_inner_radius = current_row['DN']/1000/2
        pipe_insulation_thick = current_row['insulation_thickness']
        pipe_insulation_k_value = current_row['insulation_k_value']
        
        d_pipe = {"id": str(pipe_id),
                  "node1": str(pipe_start_point),
                  "node2": str(pipe_end_point),
                  "length": str(pipe_length),
                  "innerRadius": str(pipe_inner_radius),
                  "interPipeDistance": "0.5"}
        pipe_pair = SubElement(network, "PipePair", d_pipe)

        d_pipe_properties = {
            "insulationThick": str(pipe_insulation_thick),
            "insulationkValue": str(pipe_insulation_k_value),
            "buriedDepth": "1"}
        supply_pipe = SubElement(pipe_pair, "SupplyPipe", d_pipe_properties)
        return_pipe = SubElement(pipe_pair, "ReturnPipe", d_pipe_properties)


def change_boiler_to_substation(district, substations, points, design_epsilon=0.6, design_temp_difference=20):
    print("Modifying Boilers to Substations...")

    for building in district.iter("Building"):
        # Remove boiler
        heat_source = building.find("./HeatSource")
        boiler = heat_source.find("./Boiler")
        heat_source.remove(boiler)
            
        if building.attrib["Simulate"] == "true":
            # Set Substation parameters
            building_EGID = int(building.attrib["key"])
    
            substation_row = substations.loc[substations["EGID"] == building_EGID]
            node_row = points.loc[points["EGID"] == building_EGID]

            design_factor = 115/100 # Heat exchanger sizing
            design_thermal_power = substation_row["Power[W]"].iloc[0]*design_factor   
            linked_node_id = node_row["npid"].iloc[0]

            type_substation = "simple"
            
            d_substation = {"linkedNodeId": str(linked_node_id),
                            "designThermalPower": str(design_thermal_power),
                            "designTempDifference": str(design_temp_difference),
                            "designEpsilon": str(design_epsilon),
                            "type": type_substation}
            substation = SubElement(heat_source, "Substation", d_substation)
        
        else:
            # Remove heat source
            heat_source = building.find("./HeatSource")
            building.remove(heat_source)


def add_all_dhn(district, points, pipes, substations, hs_delta_p=200000, hs_p_max = 10000000, dp_type='affine', stages = None):
    # Add District heating center
    district_heating_center = add_district_heating_center(district)
    
    # Add Network (nodes and pipes)
    add_network(district_heating_center, points=points.copy(), pipes=pipes.copy())
    
    # Add Thermal station
    ts_node_id = points.loc[points['Type']=='start heating station']['npid'].iloc[0]
    
    add_thermal_station(district_heating_center, ts_node_id, delta_p=hs_delta_p, p_max = hs_p_max, 
                        dp_type=dp_type, stages=stages, temp_storage=50, c_storage=1e7)
    
    
    change_boiler_to_substation(district, substations, points)

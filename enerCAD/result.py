# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:52:36 2023

@author: olich
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


###################################################
# 
#              Search functions
#
##################################################    

def Column_Search(df, search):
    '''
    Parameters :
    df : DataFrame to filter
    search : string to search in column names
    
    Returns :
    df_filtered : DataFrame containing filtered columns '''

    df_cols = [col for col in df.columns if search in col]
    df_filtered = df[df_cols]
    return df_filtered

###################################################
# 
#              Plots functions
#
##################################################
    
def Plot_Ideal_Annual_Demand(ideal_sum_MW): #TODO
    
    # create the ideal heating needs plot
    plt.figure('Ideal heating needs')
    plt.plot(ideal_sum_MW.index, ideal_sum_MW.values)
    plt.xlabel("Heure")
    plt.ylabel("P [MW]")
    plt.title("Demande théorique de chaleur sur l'année")
    # plt.savefig("Demande_horaire_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()

def Plot_Annual_Demand_Building(data, EGID):
    df = Column_Search(data, search="Qs(Wh)")
    df = Column_Search(df, search=str(EGID))    
    # compute Q[Wh] into P[kW], assuming constant power per hour, and knowing hourly data
    df = df/1e3
    # replace negative values with 0
    df = df.clip(lower=0)    
    # create the heating needs plot
    plt.figure(f'Heating needs EGID n°{EGID}')
    plt.plot(df.index, df.values)
    plt.xlabel("Heure")
    plt.ylabel("P [kW]")
    plt.title(f"Demande de chaleur EGID n°{EGID}")
    # plt.savefig("Demande_horaire_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()

def Plot_Annual_Demand_DHN(data):
    df = Column_Search(data, search="Qs(Wh)")
    # replace negative values with 0
    df = df.clip(lower=0) 
    df_sum = df.sum(axis=1)  
    # compute Q[Wh] into P[MW], assuming constant power per hour, and knowing hourly data
    df_sum = df_sum/1e6   
    # create the heating needs plot
    plt.figure('Heating needs')
    plt.plot(df_sum.index, df_sum.values)
    plt.xlabel("Heure")
    plt.ylabel("P [MW]")
    plt.title("Demande de chaleur du quartier")
    # plt.savefig("Demande_horaire_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()

def Plot_Internal_Gains_Building(data, EGID):
    df = Column_Search(data, search="Qi(Wh)")
    df = Column_Search(df, search=str(EGID))
    # compute Q[Wh] into P[kW], assuming constant power per hour, and knowing hourly data
    df = df/1e3
    # create the internal gains plot
    plt.figure('Internal gains')
    plt.plot(df.index, df.values)
    plt.xlabel("Heure")
    plt.ylabel("P [kW]")
    plt.title(f"Gains internes et solairesEGID n°{EGID}")
    # plt.savefig("Demande_horaire_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()

def Plot_Machine_Power_Building(data, EGID):
    df = Column_Search(data, search="MachinePower(W)")
    df = Column_Search(df, search=str(EGID))
    # compute Q[Wh] into P[kW], assuming constant power per hour, and knowing hourly data
    df = df/1e3
    # create the machine power plot
    plt.figure('Machine power')
    plt.plot(df.index, df.values)
    plt.xlabel("Heure")
    plt.ylabel("P [kW]")
    plt.title(f"Puissance des machines EGID n°{EGID}")
    # plt.savefig("Demande_horaire_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()
    

#TODO
def Plot_Machine_Power_Building_days(data, EGID):
    df = Column_Search(data, search="MachinePower(W)")
    df = Column_Search(df, search=str(EGID))    
    # compute Q[Wh] into P[kW], assuming constant power per hour, and knowing hourly data
    df = df/1e3
    # replace negative values with 0
    df = df.clip(lower=0)

    # define the four specific days
    specific_days = [79,171,264,355]
    specific_days_name = ['21 mars','21 juin','21 septembre', '21 décembre']
    
    # create the subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.flatten()
    plt.title(f"Puissance des machines EGID n°{EGID}")
    # iterate over the specific days and plot the data for 24 hours in each subplot
    for i, day in enumerate(specific_days):
        start_index = day * 24
        end_index = (day + 1) * 24    
        # extract the heating values for the specific day
        df_day = df.iloc[start_index:end_index]    
        # plot the heating values in the subplot
        axes[i].bar(df_day.index, df_day.values.flatten())
        axes[i].set_xlabel("Hour")
        axes[i].set_ylabel("kW")
        axes[i].set_title(f"Heating Needs - {specific_days_name[i]}")
        x_ticks = range(start_index+1, end_index+1)
        axes[i].set_xticks(x_ticks)
        axes[i].set_xticklabels(range(24)) 
        
    # adjust the spacing between subplots
    plt.tight_layout()
    plt.show()

def Plot_Machine_Power_DHN(data):
    df = Column_Search(data, search="MachinePower(W)")
    # replace negative values with 0
    df = df.clip(lower=0) 
    df_sum = df.sum(axis=1)  
    # compute Q[Wh] into P[MW], assuming constant power per hour, and knowing hourly data
    df_sum = df_sum/1e6 
    # create the machine power plot
    plt.figure('Machine power')
    plt.plot(df.index, df.values)
    plt.xlabel("Heure")
    plt.ylabel("P [MW]")
    plt.title("Puissance des machines du quartier")
    # plt.savefig("Demande_horaire_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()   
    
def Plot_PV_Building(data, EGID):
    df = Column_Search(data, search="SolarPVProduction(kWh)")
    # compute Q[Wh] into P[MW], assuming constant power per hour, and knowing hourly data
    df = df/1e6
    # create the machine power plot
    plt.figure('Machine power')
    plt.plot(df.index, df.values)
    plt.xlabel("Heure")
    plt.ylabel("P [MW]")
    plt.title(f"Puissance de production PV EGID n°{EGID}")
    # plt.savefig("Demande_horaire_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()

def Plot_Monotone_Demand(data):
    df = Column_Search(data, search="Qs(Wh)")
    # replace negative values with 0
    df = df.clip(lower=0) 
    df_sum = df.sum(axis=1)  
    # compute Q[Wh] into P[MW], assuming constant power per hour, and knowing hourly data
    df_sum = df_sum/1e6 
    # create the monotone de charge plot
    plt.figure('Monotone demande')
    plt.plot(range(len(df_sum)), sorted(df_sum, reverse=True))
    plt.xlabel("Heures/année")  
    plt.ylabel("P [MW]")
    plt.title("Monotone de charge de la demande")
    # plt.savefig("Monotone_de_charge_Bex_0.png", dpi=300) # save the plot as PNG with 300 dpi resolution
    plt.show()
    
def plot_load_curve(data, directory_path):
    df = Column_Search(data, search="MachinePower(W)")
    # replace negative values with 0
    df = df.clip(lower=0) 
    df_sum = df.sum(axis=1)  
    
    # create load curve plot
    plt.figure('Load duration curve')
    
    # compute P[W] into P[MW/kW], assuming constant power per hour, and knowing hourly data
    if df_sum.max()>1e6:
        df_sum = df_sum/1e6
        plt.ylabel("Power [MW]")
    else:
        if df_sum.max()>1e3:
            df_sum = df_sum/1e3
            plt.ylabel("Power [kW]")
        else:
            plt.ylabel("Power [W]")
    
    plt.plot(range(len(df_sum)), sorted(df_sum, reverse=True))
    plt.xlabel("Operating hours [h]")  
    plt.title("Buildings thermal load duration curve")
       
    load_curve_path = os.path.join(directory_path, 'Thermal_load_duration_curve.png')     
    plt.savefig(load_curve_path, dpi=300, bbox_inches='tight')
    plt.show()
    
def plot_energy_data(data, scenario_id, scenarios, directory_path):
    
    for i in range(len(scenarios)):
        scenario = scenarios[i]
        sc_id = scenario[0]
        
        if sc_id == scenario_id:
            
            production_types = scenario[2] 
            production_names = []
            colors = []
            
            for j in range(len(production_types)):
                
                prod_type = production_types[j]
                
                if prod_type == 'CHP':
                    production_names.append('Wood cogeneration')
                    colors.append('firebrick')
                if prod_type == 'Wood_boiler':
                    production_names.append('Wood boiler') 
                    colors.append('darkgoldenrod')
                if prod_type == 'Gas_boiler':
                    production_names.append('Gas boiler')
                    colors.append('black')
                if prod_type == 'Heat_Pump_Water':  
                    production_names.append('Water-to-water Heat Pump')
                    colors.append('darkturquoise')

            # Exctract the thermal station columns
            df_thermal_station = Column_Search(data, search="ThermalStation")    

            # Exctract the thermal heat production
            df_th_prod = Column_Search(df_thermal_station, search="MachinePower[") #W => Wh
            if df_th_prod.size == 0: 
                # Case where there is only one producer
                df_th_prod = Column_Search(df_thermal_station, search="MachinePower") #W => Wh
            E_max = df_th_prod.max().max()
             
            # Setup energy magnitude
            if E_max > 1e6:
                df_th_prod = df_th_prod/1e6 #MWh
                Energy_magn = 'MWh'
            else:
                if E_max > 1e3:
                    df_th_prod = df_th_prod/1e3 #kWh
                    Energy_magn = 'kWh'
                else:
                    Energy_magn = 'Wh'
                
            ### Create the stacked bar plot
            
            plt.figure(num=f'Annual Energy production scenario {scenario_id}', figsize=(10, 6),)
            
            # Plot stacked data as bars
            stacked_data = []
            for column in df_th_prod.columns:
                stacked_data.append(df_th_prod[column])

            for i, data_series in enumerate(stacked_data):
                plt.bar(df_th_prod.index, data_series, label=production_names[i], color=colors[i], alpha=1, width=1, bottom=df_th_prod.iloc[:, :i].sum(axis=1))
                        
            # Enable auto-ticker for zoomed-in x-axis
            plt.autoscale(enable=True, axis='x', tight=True)
            
            plt.xlabel('Time [hours]')
            plt.ylabel(f'Hourly satisfied thermal energy demand [{Energy_magn}]')
            plt.title(f'Annual Energy conversion to heat : scenario {scenario_id}')
            plt.legend(loc='upper center')
            plt.tight_layout()
        
            results_energy_path = os.path.join(directory_path, 'Energy_production.png')     
            plt.savefig(results_energy_path, dpi=300, bbox_inches='tight')

            plt.show()

            ### Create the operation hours pie
            
            # Count the operating status
            df_operating_hours = df_th_prod.copy()
            power_threshold = 5/1000*df_operating_hours[df_operating_hours.columns[0]].max() #threshold at 5 per 1000 of nominal power base producer
            count = 0
            for column in df_operating_hours.columns:
                if count == 0:
                    df_operating_hours[column] = df_operating_hours[column].apply(lambda x: 1 if x > power_threshold else 0)
                else:
                    df_operating_hours[column] = df_operating_hours[column].apply(lambda x: 1 if x > 0 else 0)
                count = count+1
                
            # Add column for hours without production (when base producer is off)
            df_operating_hours['Residual production'] = df_operating_hours[df_operating_hours.columns[0]].apply(lambda x: 1 if x == 0 else 0)

            # Calculate hours_on and hours_off for each machine
            machine_hours = []
            count = 0
            machine_combination = ""
            for column in df_operating_hours.columns:
                if count == 0:
                    machine_combination = production_names[count]
                    hours_on = df_operating_hours[column].sum()
                    machine = [machine_combination, hours_on]
                    machine_hours.append(machine)
                elif count == len(df_operating_hours.columns)-1:
                    hours_on = df_operating_hours[column].sum()
                    machine = ['Residual production', hours_on]
                    machine_hours.append(machine)
                else:
                    machine_combination = machine_combination+'\n+ '+production_names[count]
                    hours_on = df_operating_hours[column].sum()
                    machine = [machine_combination, hours_on]
                    machine_hours.append(machine)
                count = count + 1
                                                   
            hours = [hours for name, hours in machine_hours]
            machine_names = [name for name, hours in machine_hours]
            colors_hours = colors.copy()
            colors_hours.append('lightgrey')
            
            fig, ax = plt.subplots(num=f'Operating hours scenario {scenario_id}')
                   
            # Add a circle at the center to create the "donut hole"
            centre_circle = plt.Circle((0,0),0.20,fc='white')
            fig.gca().add_artist(centre_circle)
  
            # wedges, texts, autotexts = ax.pie(hours, labels=machine_names, colors=colors_hours, startangle=90,
            #            autopct='', pctdistance=0.85, wedgeprops=dict(width=0.3, edgecolor='w'))
            wedges, _ = ax.pie(hours, labels=None, startangle=90,
                   colors=colors_hours, wedgeprops=dict(width=0.3, edgecolor='w'))
                
            # Display actual hours inside the pie wedges
            for w, hour in zip(wedges, hours):
                ang = (w.theta2 - w.theta1) / 2. + w.theta1
                # y = 0.8 * np.sin(np.deg2rad(ang))
                # x = 0.8 * np.cos(np.deg2rad(ang))
                y = 0.9 * np.sin(np.deg2rad(ang))
                x = 0.9 * np.cos(np.deg2rad(ang))
                ax.text(x, y, f"{hour}h", ha="center", va="center", color='black', fontsize=7,
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'))

            # # Increase font size for labels and percentages
            # for text in texts + autotexts:
            #     text.set_fontsize(7)        
            
            # Equal aspect ratio ensures that pie is drawn as a circle.
            ax.axis('equal')
            
            plt.title(f"Machine Operating hours distribution : scenario {scenario_id}")
            Machine_op_hours_path = os.path.join(directory_path, 'Machine_operating_hours.png')     
            plt.savefig(Machine_op_hours_path, dpi=300, bbox_inches='tight')

            plt.show()
            

###################################################
# 
#              Exctract dataframes
#
##################################################    

def Monotone_Demand(data):
    df = Column_Search(data, search="Qs(Wh)")
    # replace negative values with 0
    df = df.clip(lower=0) 
    df_sum = df.sum(axis=1)  
    # compute Q[Wh] into P[kW], assuming constant power per hour, and knowing hourly data
    df_sum = df_sum/1e3
    df_sum = df_sum.sort_values(ascending=False).reset_index(drop=True)
    df_sum.columns = ['Number of hours [h]','Power [kW]']
    return df_sum

def get_thermal_load_curve(data):
    df = Column_Search(data, search="MachinePower(W)")
    # replace negative values with 0
    df = df.clip(lower=0) 
    df_sum = df.sum(axis=1)  
    # compute P[W] into P[kW], assuming constant power per hour, and knowing hourly data
    df_sum = df_sum/1e3
    df_sum = df_sum.sort_values(ascending=False).reset_index(drop=True)
    df_sum = df_sum.reset_index()
    df_sum.columns = ['Number of hours [h]','Power [kW]']
    return df_sum

def get_Pmax_per_EGID(data):
    df = Column_Search(data, search="MachinePower(W)")
    Power = []
    for col in df.columns:
        # Pick string ...(***)...
        EGID = col.split("(")[1].split(")")[0]
        Pmax = df[col].max()
        Power.append({"EGID":EGID,"Power[W]":Pmax})
    df_power_EGID = pd.DataFrame(Power)
    return df_power_EGID

def get_network_data_max(data):
    df = Column_Search(data, search="PipePair")
    df = Column_Search(df, search="SupplyMassFlow")
    pipe_list = []
    max_value = 0
    
    # Search for the timestep with the max value of overall network
    for column in df:
        # max_pipe = df[column].max()
        max_pipe = abs(df[column]).max()
        if max_pipe > max_value:
            max_value = max_pipe
            index_max = df[abs(df[column])==max_value].index
    index_max_value = index_max.values[0]
    
    # Extract the mass flow for the found timestep
    for column in df:
        pipe_id = column.split("PipePair")[1].split(":")[0]
        pipe_max_mass_flow = abs(df[column].loc[index_max].iloc[0])
        pipe_list.append([pipe_id,pipe_max_mass_flow])

    df_node = Column_Search(data, search="NodePair")
    df_node = Column_Search(df_node, search="SupplyTemp")
    node_list = []
    
    # Exctract the supply temperature for the found timestep
    for column in df_node:
        node_id = column.split("NodePair")[1].split(":")[0]
        node_supply_temp = abs(df_node[column].loc[index_max].iloc[0])
        node_list.append([node_id,node_supply_temp]) 
    
    df_pipe_mass_flow = pd.DataFrame(pipe_list, columns=['pid', 'Mass_flow[kg/s]']) 
    df_node_supply_temp = pd.DataFrame(node_list, columns=['npid', 'T_supply[°C]']) 
               
    return df_pipe_mass_flow, df_node_supply_temp, index_max_value

def get_energy_data(data, scenario_id, scenarios):
    
    for i in range(len(scenarios)):
        scenario = scenarios[i]
        sc_id = scenario[0]
        if sc_id == scenario_id:      
            production_types = scenario[2] 
            production_names = []
            
            for j in range(len(production_types)):
                
                prod_type = production_types[j]
                
                if prod_type == 'CHP':
                    production_names.append('Wood cogeneration')
                if prod_type == 'Wood_boiler':
                    production_names.append('Wood boiler') 
                if prod_type == 'Gas_boiler':
                    production_names.append('Gas boiler')
                if prod_type == 'Heat_Pump_Water':  
                    production_names.append('Water-to-water Heat Pump')

            fuel_consumption = []
            elec_consumption = []
            thermal_production = []
            
            ### Pump
            # Exctract the pump power consumption
            df_pump = Column_Search(data, search="PumpPower") #W
            df_pump_sum_kWh = df_pump.sum()/1e3 #kWh since hourly data
            pump_elec_consumption = df_pump_sum_kWh.values[0]
            
            ### Thermal
            df_thermal_station = Column_Search(data, search="ThermalStation")
            
            # Exctract the fuel consumption
            df_fuel = Column_Search(df_thermal_station, search="FuelConsumption[") #J
            if df_fuel.size != 0:    
                for column in df_fuel:
                    stage = column.split("[")[1].split("]")[0]
                    df_fuel_sum_kWh = df_fuel[column].sum()/3.6e6 #kWh
                    fuel_consumption.append([stage,df_fuel_sum_kWh])
            else: 
                # Case where there is only one producer
                df_fuel = Column_Search(df_thermal_station, search="FuelConsumption") #J
                df_fuel_sum_kWh = df_fuel.sum()/3.6e6 #kWh
                fuel_consumption.append([0,df_fuel_sum_kWh])
            
            # Exctract the thermal heat production
            df_th_prod = Column_Search(df_thermal_station, search="MachinePower[") #W => Wh
            if df_th_prod.size != 0:  
                # Compile yearly thermal production
                count=0
                for column in df_th_prod.columns:
                    stage = column.split("[")[1].split("]")[0]
                    df_th_prod_sum_kWh = df_th_prod[column].sum()/1e3 #kWh
                    thermal_production.append([stage,df_th_prod_sum_kWh])
                    df_th_prod = df_th_prod.rename(columns={column:production_names[count]})
                    count = count+1
            else: 
                # Case where there is only one producer
                df_th_prod = Column_Search(df_thermal_station, search="MachinePower") #W => Wh
                column = df_th_prod.columns[0]
                df_th_prod_sum_kWh = df_th_prod[column].sum()/1e3 #kWh
                thermal_production.append([0,df_th_prod_sum_kWh])
                df_th_prod = df_th_prod.rename(columns={column:production_names[0]})

            thermal_production_sum = df_th_prod.sum()/1000 #kWh   
            df_th_production_sum = pd.DataFrame(thermal_production_sum, columns=['Thermal energy production [kWh/y]'])

            ### Electrical
            # Exctract the electrical consumption
            df_el = Column_Search(df_thermal_station, search="ElectricConsumption[") #J
            if df_el.size != 0:
                for column in df_el:
                    stage = column.split("[")[1].split("]")[0]
                    df_el_sum_kWh = df_el[column].sum()/3.6e6 #kWh                
                    elec_consumption.append([stage,df_el_sum_kWh])      
                    df_el = df_el.rename(columns={column:production_names[int(stage)]})   
            else:
                # Case where there is only one producer
                df_el = Column_Search(df_thermal_station, search="ElectricConsumption") #J
                column = df_el.columns[0]
                df_el_sum_kWh = df_el[column].sum()/3.6e6 #kWh                
                elec_consumption.append([0,df_el_sum_kWh])      
                df_el = df_el.rename(columns={column:production_names[0]})
                
            electrical_consumption_sum = df_el.sum()/3.6e6 #kWh   
            df_elec_sum = pd.DataFrame(electrical_consumption_sum, columns=['Electrical energy consumption [kWh/y]'])
            
            
    return pump_elec_consumption, fuel_consumption, elec_consumption, thermal_production, df_th_production_sum, df_elec_sum




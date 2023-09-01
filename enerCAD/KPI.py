# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 17:33:39 2023

@author: Olivier Chavanne
"""

import math
import pandas as pd
import os
import matplotlib.pyplot as plt


def calculate_KPI(scenario_id, scenarios, vol_storage, df_parameters,
                  pump_consumption, fuel_consumption, elec_consumption, thermal_production):
        
    #Collect efficiencies data
    eff_columns = ['T_eff','efficiency']
    KPI_eff = df_parameters[eff_columns]
    
    #Collect technologies data
    tech_columns = ['Technology','surf_machine[m2/kW]','surf_installation[m2/kW]','lifetime[y]','cost_machine[CHF/kW]','env_machine[kgCO2/kWh]']
    KPI_tech = df_parameters[tech_columns]
    
    #Collect resources data
    resources_columns = ['Resource','cost_fuel[CHF/kWh]','env_fuel[kgCO2/kWh]']
    KPI_res = df_parameters[resources_columns]
    
    for i in range(len(scenarios)):
        scenario = scenarios[i]
        sc_id = scenario[0]
        
        if sc_id == scenario_id:
            KPI_list = []
            
            # Consumption data
            pump_elec_consumption = pump_consumption #kWh
            
            total_surface = 0
            total_cost = 0
            total_env_impact = 0
            
            # Constant values
            c_civil_eng = 300 #CHF/m3
            c_water_tank = 2000 #CHF/m3
            IRR = 5/100 # intrest rate
            O_M_rate = 5/100 # fraction of investment
            water_tank_lifetime = 25 #y
            ring_water_tank = 2 #m for maintenance
            infrastructure_h = 8 #m height
            wood_storage_h = 5 #m height
            wood_storage_inertia = 24*4 #h (4 days)
            wood_energy_volume = 4.8*650*0.65 #PCI kWh/kg *rho kg/m3 *filling
            
            # Specific values independent of technology
            Elec_buy = KPI_res[KPI_res['Resource']=='Electricity_buy'].iloc[0]
            Elec_sell = KPI_res[KPI_res['Resource']=='Electricity_CHP_sell'].iloc[0]

            c_elec_buy = Elec_buy['cost_fuel[CHF/kWh]']
            c_elec_sell = Elec_sell['cost_fuel[CHF/kWh]']

            e_elec_buy = Elec_buy['env_fuel[kgCO2/kWh]']
            e_elec_sell = Elec_sell['env_fuel[kgCO2/kWh]']
              
            ### Impact from joint facilities  
            # Water storage tank
            D_water_tank = (1.6/math.pi*vol_storage)**(1/3)
            surf_water_tank = math.pi*(D_water_tank/2+ring_water_tank)**2
            
            cost_water_tank = c_water_tank*vol_storage
            annualised_cost_water_tank = cost_water_tank*IRR/(1-1/((1+IRR)**(water_tank_lifetime)))
      
            ### Impact specific to each technology
            for j in range(len(scenario[1])):  
                unit_power = scenario[1][j]/1000 #kW
                rel_power = unit_power/sum(scenario[1])
                unit_type = scenario[2][j]
                
                fuel_cons = fuel_consumption[j][1] #kWh
                elec_cons = elec_consumption[j][1]+pump_elec_consumption*rel_power #kWh
                th_prod = thermal_production[j][1] #kWh
                
                if elec_cons > 0:
                    c_elec = c_elec_buy
                    e_elec = e_elec_buy
                else:
                    c_elec = c_elec_sell
                    e_elec = abs(e_elec_sell-e_elec_buy)
                
                if unit_type == 'CHP':
                    efficiency = 'CHP_th'
                    technology = 'CHP'
                    resource = 'Wood'
                    wood = True
                    
                if unit_type == 'Wood_boiler':
                    efficiency = 'Wood_boiler' 
                    technology = 'Wood_boiler'
                    resource = 'Wood'                    
                    wood = True
                
                if unit_type == 'Gas_boiler':
                    technology = 'Gas_boiler'
                    resource = 'Gas'                    
                    wood = False
                    
                if unit_type == 'Heat_Pump_Water':
                    technology = 'Heat_Pump_Water'
                    resource = None                
                    wood = False
                
                if wood == True:    
                    Eff = KPI_eff[KPI_eff['T_eff']==efficiency].iloc[0]
                    eff = Eff['efficiency']
                    
                Tech = KPI_tech[KPI_tech['Technology']==technology].iloc[0]
                if resource != None:
                    Res = KPI_res[KPI_res['Resource']==resource].iloc[0]
            
                ### Land use                
                # Production facilities
                s_power = Tech['surf_machine[m2/kW]']
                s_install = Tech['surf_installation[m2/kW]']
                surf_power = s_power*unit_power
                surf_hyd_el = s_install*unit_power
                surf_DHN = surf_power+surf_hyd_el
                
                # Wood storage facilities
                if wood == True:
                    volume_wood_storage = wood_storage_inertia*unit_power/eff/wood_energy_volume
                    surf_wood_storage = volume_wood_storage/wood_storage_h
                else:
                    surf_wood_storage = 0

                # Water storage tank
                rel_surf_water_tank = rel_power*surf_water_tank

                total_surf = surf_DHN+surf_wood_storage+rel_surf_water_tank
                total_surface = total_surface+total_surf

                ### Costs
                # Investment
                c_machine = Tech['cost_machine[CHF/kW]']                 
                cost_machine = c_machine*unit_power
                civil_eng_volume = surf_DHN*infrastructure_h+surf_wood_storage*wood_storage_h
                cost_civil_eng = c_civil_eng*civil_eng_volume
                investment = cost_machine+cost_civil_eng
                
                # Yearly costs
                lifetime = Tech['lifetime[y]'] 
                annuality = investment*IRR/(1-1/((1+IRR)**(lifetime)))
                cost_O_M = O_M_rate*investment
                rel_annualised_cost_water_tank = rel_power*annualised_cost_water_tank
                if resource != None:
                    c_fuel = Res['cost_fuel[CHF/kWh]']                   
                    cost_fuel = c_fuel*fuel_cons
                else:
                    cost_fuel = 0
                cost_elec = c_elec*elec_cons
                
                annual_cost = annuality+cost_O_M+rel_annualised_cost_water_tank+cost_fuel+cost_elec
                total_cost = total_cost+annual_cost
                
                ### Environmental               
                # Machines manufacturing and fuel consumption
                e_machine = Tech['env_machine[kgCO2/kWh]']
                env_machine = e_machine*th_prod/1000 #tCO2
                env_elec = e_elec*elec_cons
                
                total_env = env_machine+env_elec
                total_env_impact = total_env_impact+total_env
                
                KPI_list.append([unit_type,math.ceil(total_surf),math.ceil(annual_cost),math.ceil(total_env)])
            
            # Overall impacts
            KPI_list.append(['Total',math.ceil(total_surface),math.ceil(total_cost),math.ceil(total_env_impact)])
            
            # Create DataFrames
            df_KPI = pd.DataFrame(KPI_list, columns=['Source','land_use[m2]','cost[CHF/y]','global_warming_potential[tCO2/y]']) 

    return df_KPI

def plot_KPI(KPI_result_list, scenarios_list, directory_path):

    cost_list = []
    env_list = []
    surf_list = []
    surf_max = 0
    env_max = 0
    marker_colors = []
    marker_labels = []

    for i in range(len(scenarios_list)):
        df_KPI = KPI_result_list[i]
        
        KPI = df_KPI[df_KPI['Source']=='Total'].iloc[0]

        cost = KPI['cost[CHF/y]']
        env = KPI['global_warming_potential[tCO2/y]']
        surf = KPI['land_use[m2]']
        
        cost_list.append(cost)
        env_list.append(env)
        if env > env_max:
            env_max = env
        surf_list.append(surf)
        if surf > surf_max:
            surf_max = surf
        
        scenario = scenarios_list[i]
        if scenario == 1:
            marker_colors.append('firebrick')
        if scenario == 2:
            marker_colors.append('sienna')
        if scenario == 3:
            marker_colors.append('deepskyblue')
        if scenario == 4:  
            marker_colors.append('mediumblue')
        if scenario == 5:  
            marker_colors.append('lime')
        if scenario == 6:  
            marker_colors.append('yellowgreen')
        marker_labels.append(f'Scenario {scenario}')
    
    marker_size_max = 200
    
    x_values = cost_list
    if env_max > 1e6:
        y_values = [math.ceil(env/1e6) for env in env_list]
        env_magn = 'Mt'
    elif env_max > 1e3:
        y_values = [math.ceil(env/1e3) for env in env_list]
        env_magn = 'kt'
    else:
        y_values = env_list
        env_magn = 't'
    
    marker_sizes = [math.ceil(surf/surf_max*marker_size_max) for surf in surf_list]

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6), num='KPI Scores')
    
    # Custom graph
    ax.set_xlabel('Cost [CHF/y]')
    ax.set_ylabel(f'Global warming potential [{env_magn}CO$_{2}$eq/y]')
    ax.set_title('Key Performance Indicator Scores')
    ax.grid(True, zorder=1)
    
    # Create the scatter plot
    # scatter = ax.scatter(x_values, y_values, s=marker_sizes, alpha=1, c=marker_colors, label=f'Land use : compared to {surf_max} m$^{2}$', zorder=2)
    scatter = ax.scatter(x_values, y_values, s=marker_sizes, alpha=1, c=marker_colors, zorder=2)
        
    # Set x-axis and y-axis limits with a 20% margin
    x_margin = 0.2 * (max(x_values) - min(x_values))
    y_margin = 0.2 * (max(y_values) - min(y_values))
    ax.set_xlim(min(x_values) - x_margin, max(x_values) + x_margin)
    ax.set_ylim(min(y_values) - y_margin, max(y_values) + y_margin)
    
    # Add labels for each marker
    for i, label in enumerate(marker_labels):
        ax.text(x_values[i], y_values[i] + y_margin/4, label, fontsize=10, ha='center', va='center', color='black', zorder=3)
    # ax.legend()
    
    # Show and save the plot
    KPI_comparison_path = os.path.join(directory_path, 'KPI_comparison.png')     
    plt.savefig(KPI_comparison_path, dpi=300, bbox_inches='tight')
    plt.show()    
    

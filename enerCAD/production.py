# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:18:53 2023

@author: Olivier Chavanne
"""
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import matplotlib.pyplot as plt

###################################################
# 
#      Scenarios generation
#
##################################################

def interpolate_time(load_curve, n_hours):
    # Create a linear interpolation function
    interp_curve = interp1d(load_curve.index, load_curve["Power [kW]"], kind='linear')
    
    # Find the linear interpolation of the Power for a given Number of hours
    interpolate = int(interp_curve(n_hours))

    return interpolate

def get_scenarios(load_curve):
    
    # Power targets in [W]
    P_max = int(load_curve["Power [kW]"].max())*1000
    
    P_base = interpolate_time(load_curve, 5000)*1000
    P_base_int = int(P_max*0.8)
    P_inter = int(P_max*0.8)-P_base
    P_int_peak = P_max-P_base
    P_peak = P_max-int(P_max*0.8)
    
    # Scenario_id, Stage powers [W], Production types 
    Scenarios = [
        [1, [P_base, P_int_peak], ['CHP', 'Wood_boiler']],
        [2, [P_base, P_inter, P_peak], ['CHP', 'Wood_boiler', 'Gas_boiler']],
        [3, [P_max], ['Heat_Pump_Water']],
        [4, [P_base_int, P_peak], ['Heat_Pump_Water', 'Gas_boiler']],
        [5, [P_base, P_int_peak], ['CHP', 'Heat_Pump_Water']],
        [6, [P_base, P_inter, P_peak], ['CHP', 'Heat_Pump_Water', 'Gas_boiler']],
        [7, [P_max], ['Wood_boiler']]
        ]
    
    return Scenarios

###################################################
# 
#      Dimensionning functions
#
##################################################

def get_storage(load_curve):
    # Size the water tank storage for the Thermal station
    P_max = int(load_curve["Power [kW]"].max())*1000 #W
    E_stored = P_max*3*3600 #J, 1/8 day of storage
    
    dT = 15 #K
    cp = 4184 #J/kg/K
    rho = 1000 #kg/m3
    
    volume_storage = E_stored/(rho*cp*dT) #m3
    capacity_storage = E_stored/dT #J/K
    
    return volume_storage, capacity_storage

###################################################
# 
#      Plot result 
#
##################################################

def Plot_Scenario(df_monotone, P_dim, n_scenario):
    plt.figure()
    df_monotone["P[kW]"].plot()
    P_inf = 0
    col = ['black','blue','orange','red']
    for i in range(len(P_dim)):
        name_plot = P_dim[i][0]
        t_plot = P_dim[i][2]
        P_plot = P_dim[i][3]
        plt.plot(t_plot, P_plot, marker="o", markersize=4, markeredgecolor="red", markerfacecolor="red",label=name_plot)
        plt.plot([0, t_plot], [P_plot, P_plot], "r")
        plt.plot([t_plot, t_plot], [P_inf, P_plot], "r")
        P_inf = P_plot
    
    plt.title(f"Dimensionnement de la production : scénario {n_scenario}")
    plt.xlabel("Heures classées")
    plt.ylabel("P [kW]")
    plt.legend(loc="upper right")
    




    
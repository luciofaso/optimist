# -*- coding: utf-8 -*-
"""
    Created on Fri Jun  9 14:40:46 2017
    @author: lraso

"""

## TEST
import optimist
import pandas as pd

#%% import inflow data

# Load CSV file
Location='../Data/'
csv_file_name_1='HydrologicalRawData.csv'
discharge=pd.read_csv(Location+csv_file_name_1,sep=',',dayfirst=True, parse_dates=['Date'],index_col='Date')


#%% create system
system_0=System('test')


#%% create basin (inflow to the reservoir)
bafing=Inflow('Bafing',senegal_river,discharge['Bafing'])

#%% create simplified reservoir
V_max=10^13

res=Reservoir('test res',system_0,V_max,releases='Turbine',evaporation=0,losses=20,inflow=bafing)#


#OR res.add_input(bafing)


def rule_1(vol):
    release=0.01*vol/(24*60)
    return(release)

res.release_rules=rule_1


#%% SIMULATION

#initialize all variables
for element in system_0.internal_elements:
    for var in element.variables:
        var['init']= 0

# or res.set_initial_condition(V_max) 

system_0.sim_H()
 

#%% check results


  
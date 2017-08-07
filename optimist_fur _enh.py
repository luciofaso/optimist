# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:46:30 2017

@author: lraso and B.Dieu
"""
from pympler.tracker import SummaryTracker
#tracker = SummaryTracker()
#tracker.print_diff()



global count 
count = 0

#%% importing required modules 
import numpy as np #use log, matrix, matrix operations, (multivariate) normal ditribution
import pandas as pd
#import some_optimization_package as short_name #variable and constraints

from scipy.interpolate import interp1d

#part optimize

from scipy.optimize import brentq
import math
import matplotlib.pyplot as plt

#part platypus

from platypus import NSGAII,Problem,Real,Evaluator, unique, nondominated
####
import matplotlib.pyplot as plt
#####

import os
import time

# Small change to test github

# %% class Element


class Element:
    """ 
    
    Abstract class Element
    
    _Element is the basic class, xxx for both `Element` and `System`
    
            
        
    Attributes
    ----------
    :system: Element
        is the system to  which the element belongs
        NB : itself in case of the system
        
    :internal_elements: list
        is a list of the elements of the system 
        
    :variables: list
        is a list of global (for the system) variables
        
    :name: str
        is the name of the element 
    
    """
    
    def __init__(self,name,system = 'itself',variables=[],**kwargs): 
        """
        
        Create _Element object
    
        :system: Element
                The system to which the element belongs (same class)
                
        :variables: list
                    list of variable's names
                    
        :kwargs: is unused
        
        """ 
        if system == 'itself': #case where we initialize the system
            self.system = self
            self.internal_elements = []
            self.variables = variables
        else:                  #case where we add an internal element
            self.system = system 
            self.system.internal_elements += [self]
            
        self.name = name         
        
    def __repr__(self):
        """
        
        Display the class, the name , and all the names of the object's 
        variable
        
        """
        print('Class: ' + self.__class__.__name__)
        print('Name: ' + self.name+'\n')
        print('List of variables')
        for var in self.variables :
            print(var['name']) 

    def add_variable(self,name,min_value=None,max_value=None,initial_value=None):
        """
        
        This method add a variable to an object.
        
        :self: Element
            is the object that owns the variable
            
        :name: str
            is the name of the variable 
            
        :min_value: float
            is the minimal value of the variable
              NB: this feature isn't exploited yet but should be in a wieder 
                  application
                  
        :max_value: float
            is the minimal value of the variable
              NB: this feature isn't exploited yet but should be in a wieder 
                  application
                  
        :initial_value: float
            is the initial value of the variable
        
        
        STRUCTURE
        ---------
        
        The structure of every variable is the following. The "variables" is a
        dict which contains keys corresponding to :
        :'name': str
             the name of the variable 
             
        :'min_value': float 
            is the minimal value of the variable
              NB: this feature isn't exploited yet but should be in a wieder 
                  application
                  
        :'max_value': float
            is the minimal value of the variable
              NB: this feature isn't exploited yet but should be in a wieder 
                  application
                  
        :values: list
            is the list of the values of the variable
        
'        """
        # Test to ensure the unicity of the name of  the variable at the 
        # system's level and at the object level as well.
        var_in_sys = any (var['name'] == name for var in self.system.variables)
        var_in_obj = any (var['name'] == name for var in self.variables)  
        # NB : Test above done two times  if  self==system
        if var_in_sys or var_in_obj:
            print('This variable already exists')
        else:
            variable = {
                 'name':name,
                 'min_value':min_value,
                 'max_value':max_value,
                 'values':[None]
            }
            variable['values'][0] = initial_value 
            self.variables +=  [variable] 

    def set_initial_conditions(self,name_variable,initial_value):
        """
        
        Set initial condition (initial value) for variable to initial_value
        
        :self: Element
            is the object that owns the variable
            
        :name_variable: str
            is the name of the variable
            
        :initial_value: float
            is the value used for the initialization
                
        """
        
        missing = True
        for var in self.variables:
            if var['name']==name_variable:
                missing = False
#                if var['values'][0] != None:
#                    print('!WARNING! '+name_variable+' already initialized')
                var['values'][0] = initial_value 
        if missing:
            print('!ERROR! '+name_variable+' does  not exist')
            
    def write_variable(self,name_variable,t,value):
        """
        
        This method modifies a specific value of a variable
        
        :self: Element
            is the object which owns the variable
            
        :name_variable: str
            is the name of the variable
            
        :t: int
            is the index of the value
            
        :value: float
            is the value to write
        
        """
        missing = True
        for var in self.variables:
            if var['name'] == name_variable:
                missing = False
                var['values'][t]= value
        if missing:
             print('!!ERROR!! '+name_variable+' does  not exist')
             
    def init_sim(self,df_input):
        """
        
        This method initializes the simulation, which means it sets the lenght
        of the values of each variables in order to match with the duration of
        the simulation.
        It also sets by computation the first value for some variables.
        
        :self: Element
            has to be a system
            
        :df_input: array
            is an array of values which conditiones the length of the 
            simulation.
        
        """
        
        H = len(df_input[0])
            #create empty lists for all variables:
        for var in self.variables:
            var['values'][1:H] = [None] * (H-1)
        for  elt  in self.internal_elements: 
            if type(elt) == Inflow:
                pass
            else:    
                for var  in  elt.variables:
                    if var['name'] == 'water_altitude_after':
                        elt.set_initial_conditions(
                            var['name'],
                            altitude_after_manantali(0)
                            )
                    if var['name'] == 'water_altitude_before':
                        vol = get_value(elt,'volume',0)
                        elt.set_initial_conditions(
                            var['name'],
                            altitude_before_manantali(vol)
                            )
                    var['values'][1:H] = [None] * (H-1)
   
    def _sim_step(self,t,kwargs):
        """ 
        Run one step of simulation for all internal elements 
        
        :self: Element
            !!! has to be a system !!!
            
        :t: int
            is the index of the step to be simulated
            
        :kwargs: dict
            contains the parameters for the rule
        
        """
        for element in self.internal_elements:
            element._sim(t,kwargs)
            
        
        
    def sim_H(self,kwargs):
        """
        
        Compute the simulation for the horizon H, for each time t calls _sim
        
        :self: Element
            !!! has to be a system !!!
            
        :kwargs: dict
            contains the parameters for the rule
        
        """        
        b.set_initial_conditions('volume',b.max_volume)
        df_input = input_test 
        
        H = len (df_input[0])
        
        self.init_sim(df_input)
            
        for t in range(1,H):
            self._sim_step(t,kwargs)
            
        
class Reservoir(Element):
    """ 
    
    Reservoir element
    
    Attributes
    ----------    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
        
    :name: str
        is the name of the element 
    
    :name_upper_elt: str
        is the name of the element situated directly upstream
        It defines the inflow received by the reservoir
        
    :max_volume: float
        Reservoir maximum volume [:math:`m^3`]
    
            
            
    Variables
    ---------
    
    :volume:
        is the volume of water in the reservoir
        
    :spillage:
        is the discharge of spilled water 
        
    :turbined_release:
        is the discharge of turbined water 
        
    :outflow:
        is the total discharge of water 
        
    :water_altitude_before:
        is the altitude of the water in the reservoir
        
    :water_altitude_after:
        is the altitude of the water on the river right after the reservoir
        
    :produced_energy:
        is the amount of energy produced    
    
    """
    def __init__(self,
            name,
            system,
            name_upper_elt,
            max_volume,
            max_turbined,
            max_spillage,
            **kwargs
            ):
        """ """
        self.system = system
        self.name = name
        self.name_upper_elt = name_upper_elt
        self.system.internal_elements += [self]
        self.variables = []
        self.max_volume = max_volume
        min_volume= kwargs['min_volume'] if 'min_volume' in kwargs else 0
        self.coeff_prod = kwargs['coeff_prod'] if 'coeff_prod' in kwargs else 1
        self.add_variable('volume',min_volume,max_volume)
        self.max_spillage = max_spillage
        self.add_variable('spillage',0,max_spillage,0)
        self.max_turbined = max_turbined
        self.add_variable('turbined_release',0,max_turbined,0)
        self.add_variable('water_altitude_before',0,None,205)
        self.add_variable('water_altitude_after',0,None,153)
        self.add_variable('produced_energy',0,None,0)
        self.add_variable('outflow',0,None,0)
          
    def _sim(self,t,kwargs):
        """
        
        Compute one step of the simulation

        :t: int
            is the index of the step to be simulated
            
        :kwargs: dict
            contains the parameters for the rule
        
        """
        x_t_min_1 = get_value(self,'volume',t-1)        
        inflow = get_value(
                     get_component(self.system,self.name_upper_elt),
                     'outflow',
                     t
                     )
        releases=RR(self,t,kwargs)
                
        evap = tabl_evap[t%365]
        turb = releases[0]
        spill = releases[1]
        self.write_variable('turbined_release',t,turb)
        self.write_variable('spillage',t,spill)
        outflow = turb+spill-evap

        self.write_variable('outflow',t,outflow)

        new_volume = x_t_min_1+(86400*(-outflow+inflow))
        
        self.write_variable('volume',t,new_volume)
        alt_be = altitude_before_manantali(new_volume)
        alt_af = altitude_after_manantali(outflow)
        self.write_variable('water_altitude_before',t,alt_be)
        self.write_variable('water_altitude_after',t,alt_af)
        prod_energy = 1000*9.81*0.84*(24*10**-6)*(alt_be-alt_af)*turb
        self.write_variable('produced_energy',t,prod_energy)
        
                 
class Confluence(Element):
    """ 
    
    The purpose of this object is to add the outflow of different elements to 
    merge them into one unique outflow.
    
    Attributes
    ----------    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
        
    :name: str
        is the name of the element 
    
    :list_name_upper_elt: list of str
        is a list of the names of the element situated directly upstream
        It defines the inflow received by the confluent
        
    :delay: list of float
        is the list containing the times for the water to flow from each 
        upstream element to the point of confluence
        
    :coeff: list of float
        is the list containing coefficients which traduce the evolution of the
        discharge considering hydrologic and topologic phenomena.
     
    Variables
    ---------
            
    :outflow:
        is the total discharge of water 
        
    """
    def  __init__(self,name,system,list_name_upper_elts):
        """ 
        
        Initialization of the object
        
        :system: Element
            The system to which the Reservoir belongs to
            !!! has to be a system !!!
        
        :name: str
            is the name of the element 
    
        :list_name_upper_elt: list of str
            is a list of the names of the element situated directly upstream
            It defines the inflow received by the confluent
            
        """
        self.system = system
        self.name = name
        self.list_name_upper_elts = list_name_upper_elts
        self.system.internal_elements += [self]
        self.variables = []
        self.delay = [4,4,3]
        self.coeff = [1.36,1.36,1.19]
        self.add_variable('outflow',0,None,0)
        self.add_variable('St Louis outflow',0,None,0)
        self.add_variable('avg_outflow')
        
    def actualise_flow(self):
        """
        
        This method computes the sum of the inflows and stock the values in 
        the variable inflow for the whole horizon of the simulation.
        This is why it is called at the end.
        
        """
        name_elt = self.list_name_upper_elts[0]
        up = get_component(self.system,name_elt)
        inflow = get_variable(up,'outflow')['values']
        outflow = get_variable(self,'outflow')['values']
        outflow[0:len(inflow)]= [0]*len(inflow)
        max_delay = max(self.delay)
        for t in range(max_delay):
            outflow[t] = 0
        for i in range(len(self.list_name_upper_elts)):
            name_elt = self.list_name_upper_elts[i]
            delay = self.delay[i]
            coeff = self.coeff[i]
            up = get_component(self.system,name_elt)
            inflow = get_variable(up,'outflow')['values'] 
            for t in range(len(inflow)-max_delay):
                outflow[t+delay] += coeff*inflow[t]
        get_variable(self,'outflow')['values'] = outflow 

                    
    def actualise_flow_stls(self):
        """
        
        This method computes the discharge at Saint Louis taking into account
        the evaporation on the river.
        
        to be called after actualise_flow
        
        """
        outflow = get_variable(self,'outflow')['values']
        evap = [156,184,234,234,204,98,0,0,0,48,96,99]
        outflow_stls = outflow
        for t in range(len(outflow)):
            outflow_stls[t] = outflow[t]-evap[(t%365)//31]
        get_variable(self,'St Louis outflow')['values'] = outflow_stls
            
    def  create_moving_average(self):
        """
        
        This method computes the moving averages of the height of the water.
        It creates a new variable for each width and stocks the values.
        
        :windows: list of int
            it contains the different values for the width of the windows for 
            the moving average.
        It is called after actualise_height.   
        
        """
        result = moving_average(
                     get_variable(self,'outflow')['values'],
                     5
                     )
        get_variable(self,'avg_outflow')['values'] = result
                    
    
    def _sim(self,no,none):
        pass

        
class FloodPlain(Element):
    """
    This element represents a plain where the flood agriculture is used.
    Considering the inflow and other characteristics it provides us the flooded
    area for each  year.

    Attributes
    ----------    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
        
    :name: str
        is the name of the element 
    
    :name_upper_elt: str
        is the name of the element situated directly upstream
        It defines the inflow received by the plain.
        
    Variables
    ---------
            
    :height_bakel: 
        is the height of the water at Bakel
        
    :inw_avg%?%:
        
    
    """
    def  __init__(self,name,system,name_upper_elt):
        """ """
        self.system = system
        self.name = name
        self.name_upper_elt = name_upper_elt
        self.system.internal_elements += [self]
        self.variables = []
        self.add_variable('height_bakel',0)
        for  size in [6,10,15,20,25,30] :
            self.add_variable('height_avg'+str(size))
        
    def _sim(self,none1,none2):
        pass
        
    def actualise_height(self):
        """
        
        This method computes the height of the water at Bakel and stocks the 
        values in the variable height_bakel for the whole horizon of the 
        simulation.
        This is why it is called at the end.
        
        """
        up = get_component(self.system,self.name_upper_elt)
        inflow = get_variable(up,'outflow')['values']
        l = len(inflow)
        for t in range(l):
            self.write_variable('height_bakel',t,height_bakel(inflow[t]))
            
    def  create_moving_averages(self,windows=[6,10,15,20,25,30]):
        """
        
        This method computes the moving averages of the height of the water.
        It creates a new variable for each width and stocks the values.
        
        :windows: list of int
            it contains the different values for the width of the windows for 
            the moving average.
        It is called after actualise_height.   
        
        """
        for  size in windows :
            result = moving_average(
                         get_variable(self,'height_bakel')['values'],
                         size
                         )
            get_variable(self,'height_avg'+str(size))['values'] = result
        
               
class Inflow(Element):
    """
    
    Can be considered as a water spring

    Attributes
    ----------    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
        
    :name: str
        is the name of the element 
    
    :values: array-like
        is the values of the flow for the whole horizon.
    
            
            
    Variables
    ---------
        
    :outflow:
        is the total discharge of water 
    
    """
    def  __init__(self,name,system,values):
        self.system = system
        self.name = name
        self.system.internal_elements += [self]
        self.values = values
        self.variables = []
        self.add_variable('outflow',0,None,0)
        H = len(values)
        get_variable(self,'outflow')['values'][1:H] = [None] * (H-1)
        for t in range(H):
            self.write_variable('outflow',t,values[t])
            
    def _sim(self,no,none):
        pass
    

def ind_mean_deficit(system,component_name,variable_name,floor):
    """
    
    This function is used to compute the values of an indicator for each year 
    of the horizon.
    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
     
    :component_name: str
        is the name of the object which owns the variable
            
    :variable_name: str
        is the name of the variable
        
    :floor:
        is the minimal value, under which the deficit begins
    Output
    ------
    
    :indicator: list of floats
        a list containing the values of an indicator for each year 
    of the horizon.
            
    """
    data = get_variable(get_component(system,component_name),variable_name)
    data = data['values']
    l = len(data)
    nbr_year = int(np.floor(l/365))
    indicator = []
    for  i in range(nbr_year):
        mean = 0
        if type(floor) == int :
            for val in data[365*i:365*(i+1)]:
                mean +=  max(0,floor-val)
            indicator += [mean/365]
        else:
            for j in range(365*i,365*(i+1)):
                mean +=  max(0,floor[j]-data[j])
            indicator += [mean/365]
    return indicator

def ind_mean_excess(system,component_name,variable_name,ceil):
    """
    
    This function is used to compute the values of an indicator for each year 
    of the horizon.
    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
     
    :component_name: str
        is the name of the object which owns the variable
            
    :variable_name: str
        is the name of the variable
        
    :ceil: float
    is the maximum value, above which the excess begins
        
    Output
    ------
    
    :indicator: list of floats
        a list containing the values of an indicator for each year 
    of the horizon.
            
    """
    data = get_variable(get_component(system,component_name),variable_name)
    data = data['values']
    l = len(data)
    nbr_year = int(np.floor(l/365))
    indicator = []
    for  i in range(nbr_year):
        mean = 0
        for val in data[365*i:365*(i+1)]:
            mean +=  max(0,val-ceil)
        indicator += [mean/365]
    return indicator

def ind_mean(system,component_name,variable_name):
    """
    
    This function is used to compute the values of an indicator for each year 
    of the horizon.
    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
     
    :component_name: str
        is the name of the object which owns the variable
            
    :variable_name: str
        is the name of the variable
        
    Output
    ------
    
    :indicator: list of floats
        a list containing the values of an indicator for each year 
    of the horizon.
            
    """
    data = get_variable(get_component(system,component_name),variable_name)
    data = data['values']
    l = len(data)
    nbr_year = int(np.floor(l/365))
    indicator = []
    for  i in range(nbr_year):
        indicator += [np.mean(data[365*i:365*(i+1)])]
    return indicator

def ind_count_deficit(system,component_name,variable_name,minimal_value):
    """
    
    This function is used to compute the values of an indicator for each year 
    of the horizon.
    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
     
    :component_name: str
        is the name of the object which owns the variable
            
    :variable_name: str
        is the name of the variable
        
    :minimal_value: float
        is the value under which we consider the condition violated 
        
    Output
    ------
    
    :indicator: list of floats
        a list containing the values of an indicator for each year 
    of the horizon.
            
    """
    data = get_variable(get_component(system,component_name),variable_name)
    data = data['values']
    l = len(data)
    nbr_year = int(np.floor(l/365))
    indicator = []
    for  i in range(nbr_year):
        counter = 0
        for val in data[365*i:365*(i+1)]:
            if val < minimal_value :
                counter += 1
        indicator += [counter]
    return indicator   
    
def  moving_average(data,size):
    ''' data is a list
        size is the size of the moving window
    '''
    L = len(data)
    if size > L :
        print('!ERROR! Window too wide')
        return None
    result =[0]*L
    summ = 0             
    for i in range(0,size):
        summ = summ + data[i]
        result[i] = summ / (i+1)
    for i in range(size,L):
        summ = summ - data[i-size] + data[i]
        result[i] = summ / size
    return result

def ind_flood(system,component_name):
    """
    
    This function is used to compute the values of an indicator for each year 
    of the horizon.
    
    :system: Element
        The system to which the Reservoir belongs to
        !!! has to be a system !!!
     
    :component_name: str
        is the name of the object which owns the variable
            
    Output
    ------
    
    :indicator: list of floats
        a list containing the values of an indicator for each year 
    of the horizon.
            
    """
    find_name = []
    for  var in get_component(system,component_name).variables:
        if var['name'][0:10] == 'height_avg':
            data = var['values']        
            l = len(data)
            nbr_year = int(np.floor(l/365))
            maxes_for_each_years = []
            for  i in range(nbr_year):
                maxes_for_each_years += [height_bakel(max(data[365*i:365*(i+1)]))]
            find_name += [maxes_for_each_years]
    indicator = []
    sizes = [6,10,15,20,25,30]
    for  i in range(len(find_name[0])):
        min_area = area_flood(find_name[0][i],sizes[0])
        for j in range(6): 
            area = area_flood(find_name[j][i],sizes[j])
            if area < min_area:
                min_area = area
        indicator += [min_area]
    return indicator
       
def get_value(component,variable_name,t):
    """
    
    This function returns the value of a variable at one precise point.
    
    :component: Element
        is the object which owns the variable
    :variable_name: str
        is the name of the variable
    :t: int
        is the index of the wanted value.
        
    """
    missing = True
    for var  in component.variables:
       if var['name']==variable_name:
            missing = False
            result = var['values'][t]
    if missing:
         print('!ERROR! '+variable_name+' does  not exist')
         return None
    else:
         return result
        
def get_component(system,component_name):
    """
    
    This function returns the value of a variable at one precise point.
    
    :component: Element
        is the object which owns the variable
    :variable_name: str
        is the name of the variable
    :t: int
        is the index of the wanted value.
        
    """
    missing = True
    for elt  in system.internal_elements :
       if elt.name == component_name:
            missing = False
            result = elt
    if missing:
        print('!ERROR! '+component_name+' does  not exist')
        return None
    else:
        return result        
        
def get_variable(component,variable_name):
    """
    
    This function returns a variable.
    
    :component: Element
        is the object which owns the variable
    :variable_name: str
        is the name of the variable
        
    """
    missing = True
    for var  in component.variables:
       if var['name']==variable_name:
            missing = False
            result = var
    if missing:
        print('!ERROR! '+variable_name+' does  not exists')
        return None
    else:
        return result        
          
def RR(component,t,kwargs) : # system  e.g.the reservoir
    """
    
    Linked with an operative element (a component),
       Uses the inputs to compute and return the command.
       Inputs
    ---------
    :system: System
            The system to which the Release rule is linked
            Instance of element
            
    :in kwargs : parameters needed to define the shape and the nature of
                   the rule
                   dict containing whose values are floats 
    
    """
    releases = rule(component,t,kwargs)
    return releases
                                          
def rule(component,t,kwargs) :
    """  compute the releases"""


    inflow = input_test[0][t]
    v_t_1 = get_value(component,'volume',t-1)
    c1 = [kwargs['c1_1'],kwargs['c2_1'],kwargs['c3_1'],kwargs['c4_1']]
    c2 = [kwargs['c1_2'],kwargs['c2_2'],kwargs['c3_2'],kwargs['c4_2']]
    r = [kwargs['r1'],kwargs['r2'],kwargs['r3'],kwargs['r4']]
    w = [kwargs['w1'],kwargs['w2'],kwargs['w3'],kwargs['w4']]
    data = [v_t_1,t]
    release = max(a_t(data,c1,c2,r,w,component.max_spillage),continuity_mins[(t%365)//31])
    #case for the flood
    if altitude_before_manantali(v_t_1) > fct_limit_height(t%365) :
        turb_effect = component.max_turbined
        spill_effect = component.max_spillage
    else:
        
        turb_effect =min((v_t_1/86400)+inflow,release,component.max_turbined)
        spill_effect = min(
                           (v_t_1/86400)+inflow-turb_effect,
                           release-turb_effect,
                           component.max_spillage
                           )
    return [turb_effect,spill_effect]

######### Data extraction
def extract():
    ## Load CSV file
    Location = ''
    CsvFileName = 'HydrologicalRawData.csv'
    DischargeData = pd.read_csv(
                        Location+CsvFileName,
                        header=0,
                        sep=',',
                        dayfirst=True,
                        parse_dates=['Date'],
                        index_col='Date'
                        )
    Bafing = DischargeData['Bafing']
    Oualia = DischargeData['Oualia']
    Kidira = DischargeData['Kidira']
    Bafing_ar = Bafing.values
    Oualia_ar = Oualia.values
    Kidira_ar = Kidira.values
    print("extracted")
    return [Bafing_ar,Oualia_ar,Kidira_ar] 
    #####!!!!!!!!!just for  tests
    
def evap_interp():
    """
    The function returns the value of the evaporation for Manantali dam for 
    every day.
    
    """
    X = range(1,366)
    XP = np.linspace(15,365-15,12)
    Y = [156,184,234,234,203,98,-92,-181,-60,48,96,99]
    for val in Y:
        val = val*4*10**8/(30.4*1000*86400)
        res = np.interp(X,XP,Y,None,None,365)
    return res

def fct_limit_height(t):
    if t <= 200:
        return 210.5
    elif t <= 230:
        return -0.2*(t-200)+210.5
    elif t <= 285:
        return 0.101*(t-230)+204.5
    else:
        return 210.5
    
def creation_limit_height():
    limit_height = [0]*365
    for t in range(365):
        limit_height[t] = fct_limit_height(t)
    return limit_height
############ Physical relations
def height_bakel(q):
    """
    returns the height at bakel in cm  knowing  the discharge (m**3/s)
    :q: float
         discharge 
    
    Output
    ------
    
    height at bakel float
    
    """
    return 0.0037*(q**2)-0.0917*q+35.22
    
def area_flood(h,N):
    """
    
    Returns the flooded area in ha depending on the size of the moving
    average N and the highest value of water height
    
    :h: float
        highest value of water at bakel
    :N: int
        is the width of the window
        
    Output
    ------
    
    :area: float
        the value of the area flooded.    
    
    """
    if N==6:
        value = 212.34*h-110604
        if value >= 45000:
            return value
        else:
            return 0
    elif N==10:
        value = 212.62*h-106041
        if value >= 45000:
            return value
        else:
            return 0
    elif N==15:
        value = 221.79*h-107845
        if value >= 45000:
            return value
        else:
            return 0
    elif N==20:
        value = 228.43*h-107578
        if value >= 45000:
            return value
        else:
            return 0
    elif N==25:
        value = 236.47*h-108188
        if value >= 45000:
            return value
        else:
            return 0
    elif N==30:
        value = 240.64*h-106174
        if value >= 45000:
            return value
        else:
            return 0
    else:
        print('!!ERROR!! Invalid value for N')
        return None
      
def altitude_after_manantali(q):
    """
    Compute the altitude (in m) of the water right after Manantali dam,
    knowing the discharge (m**3/s)
    
    """
    return -3*(10**(-7)*(q**2))+0.0031*q+153.82
               
def altitude_before_manantali(v):
    """
    Compute the altitude  (in m)of the water right before  Manantali dam,
    knowing the  volume (m**3)
    
    """
    return 4*(10**-29)*(v**3)-(10**-18)*(v**2)+(10**-8)*v+157
               
############## optimisation 

def GP_sim(system,kwargs):
    """
    
    is the function used for the optimization, the one called several times by
    platypus
    
    """
    system.sim_H(kwargs)
    d.actualise_flow()
    c.actualise_height()
    c.create_moving_averages()
    d.actualise_flow()
    d.actualise_flow_stls()
    d.create_moving_average()
    
#    ##ECOLOGY
#    J1 = ind_mean_deficit(a,'Fusion','St Louis outflow',14)
#    J1M = np.mean(J1)
#    J1D = np.percentile(J1, np.arange(0, 100, 10))[-1]
    ##Flood Defense
    J2 = ind_mean_excess(a,'Fusion','outflow',4500)
    J2M = np.mean(J2)
    J2P = np.percentile(J2, np.arange(0, 100, 10))[-1]
#    ##Irigation
#    J3 = ind_mean_deficit(a,'Fusion','avg_outflow',10)
#    J3M = np.mean(J3)
#    J3P = np.percentile(J3, np.arange(0, 100, 10))[-1]
    ##NAVIGATION
    J4 = ind_count_deficit(a,'Fusion','outflow',10)
    J4M = np.mean(J4)
    J4P = np.percentile(J4, np.arange(0, 100, 10))[-1]
    ##ENERGY PRODUCTION
    J5 = ind_mean(a,'Manantali','produced_energy')
    J5M = np.mean(J5)
    J5P = np.percentile(J5, np.arange(0, 100, 10))[0]
    ##FLOOD PRODUCTION
    J6 = ind_flood(a,'Plain')
    J6M = np.mean(J6)
    J6P = np.percentile(J6, np.arange(0, 100, 10))[0]
#    print(get_value(b,'turbined_release',150))
    return([J2M,J4M,J5M,J6M])  

def a_t(X, 
        c1=[],
        c2=[],
        r=[],
        w=[],
        max_val_release=1,
        n=4
        ):
    a = 0
    for j in range(n):
        dist_1=(X[0]-c1[j])**2
        dist_2=(X[1]-c2[j])**2
        a += max_val_release*(w[j]*(math.exp(-((dist_1+dist_2)/r[j]**2))))

    return a

############## graphic
  
def plot(component,variable_name, window=None):
    """
    
    """
    values = get_variable(component,variable_name)['values']
    if window==None:
        plt.plot(values)
        plt.ylabel(variable_name+' of '+component.name)
        plt.show()
        
def plot_ind():

    X =[]
    Y =[]
    data = unique(nondominated(algorithm.result))
    for solution in data:
        X =[]
        Y =[]
        for i in range(12):
            if i == 10 or i== 11:
                X += [i]
                Y += [solution.objectives[i]/(10**7)]
            elif i == 8 or i== 9:
                X += [i]
                Y += [solution.objectives[i]/(10**2)] 
            else:
                X += [i]
                Y += [solution.objectives[i]]
        plt.plot(X,Y)
#    Coor = sorted(Coor,key=lambda l : l[0])
#    X = []
#    Y = []
#    for elt  in Coor :
#        X += [elt[0]]
#        Y += [elt[1]]

#    plt.scatter(X,Y)
    plt.plot(X,Y)
    plt.show()
   

class Lake (Problem):
    """Definition of the Lake's problem"""
    def __init__(self,nbr_eval) :
        super(Lake,self).__init__(16,4,1)
        self.types[:] = [Real(c1_1_min,c1_1_max),
                         Real(c1_2_min,c1_2_max),
                         Real(c2_1_min,c2_1_max),
                         Real(c2_2_min,c2_2_max),
                         Real(c3_1_min,c3_1_max),
                         Real(c3_2_min,c3_2_max),
                         Real(c4_1_min,c4_1_max),
                         Real(c4_2_min,c4_2_max),
                         Real(r1_min,r1_max),
                         Real(r2_min,r2_max),
                         Real(r3_min,r3_max),    
                         Real(r4_min,r4_max), 
                         Real(0,1),
                         Real(0,1),
                         Real(0,1),
                         Real(0,1)
                         ]
        #TODO no problem with w2 and w1 ????
        self.directions[:] = [Problem.MINIMIZE, 
#                              Problem.MINIMIZE, 
#                              Problem.MINIMIZE,                               
                              Problem.MINIMIZE, 
                              Problem.MAXIMIZE, 
                              Problem.MAXIMIZE]
        self.constraints[:] = "==0"
        self.count = 0  ###personal
        self.nbr_eval = nbr_eval  ###personal
    def evaluate(self,solution):
        VAR = solution.variables
        ind_test = GP_sim(a,
                          {'c1_1' : VAR[0],
                           'c1_2' : VAR[1],
                           'c2_1' : VAR[2],
                           'c2_2' : VAR[3],
                           'c3_1' : VAR[4],
                           'c3_2' : VAR[5],
                           'c4_1' : VAR[6],
                           'c4_2' : VAR[7],
                           'r1' : VAR[8],
                           'r2' : VAR[9],
                           'r3' : VAR[10],
                           'r4' : VAR[11],
                           'w1' : VAR[12],
                           'w2' : VAR[13],
                           'w3' : VAR[14],
                           'w4' : VAR[15]
                           }
                          )
        self.count += 1
        progress = ((self.count)*100)/(self.nbr_eval)
        if progress%1 == 0:
            print(str(progress)+'%')
        solution.objectives[:] =  ind_test
        solution.constraints[:]=[VAR[12]+VAR[13]+VAR[14]+VAR[15]-1]
                           
                           
########## Call of extraction
continuity_mins = [130,150,200,215,230,250,250,200,170,150,130,130]

data = extract()

nbr_year = 50
if type(nbr_year)==int:
    input_test = [data[0][0:nbr_year*365],
                  data[1][0:nbr_year*365],
                  data[2][0:nbr_year*365],
                 ]
else:
    input_test = data


limit_height = creation_limit_height()

tabl_evap = evap_interp()


######### Definition of the system

a = Element('Senegal River')
a1 = Inflow('Bafing',a,input_test[0])
a2 = Inflow('Oualia',a,input_test[1])
a3 = Inflow('Kidira',a,input_test[2])


b = Reservoir('Manantali',a,'Bafing',10**10,500,5500,coeff_prod = 20)

d = Confluence('Fusion',a,['Manantali','Oualia','Kidira'])



c = FloodPlain('Plain',a,'Fusion')

#print(GP_sim(a,{'c1' : 10000000,
#                'c2' : 10000000,
#                'c3' : 10000000,
#                'c4' : 10000000,
#                'r1' : 10000000000,
#                'r2' : 10000000000,
#                'r3' : 10000000000,
#                'r4' : 2000000000,
#                'w1' : 0.5,
#                'w2' : 0.5,
#                'w3' : 0.5,
#                'w4' : 0.5
#                }
#             ))



c1_1_min = 1
c1_1_max = 10000000000
c1_2_min = 0
c1_2_max = 365
c2_1_min = 1
c2_1_max = 10000000000
c2_2_min = 0
c2_2_max = 365
c3_1_min = 1
c3_1_max = 10000000000
c3_2_min = 0
c3_2_max = 365
c4_1_min = 1
c4_1_max = 10000000000
c4_2_min = 0
c4_2_max = 365
r1_min = 10
r1_max = 20000000000
r2_min = 10
r2_max = 20000000000
r3_min = 10
r3_max = 20000000000
r4_min = 10
r4_max = 20000000000


parameters = [c1_1_min,c1_1_max,c1_2_min,c1_2_max,
              c2_1_min,c2_1_max,c2_2_min,c2_2_max,
              c3_1_min,c3_1_max,c3_2_min,c3_2_max,
              c4_1_min,c4_1_max,c4_2_min,c4_2_max
              ]
parameters += [r1_min,r1_max,r2_min,r2_max,r3_min,r3_max,r4_min,r4_max,]


nbr_eval = 1000

########## Running optimisation

class_for_opt = Lake(nbr_eval)
algorithm = NSGAIII(class_for_opt)
start = time.clock()
algorithm.run(nbr_eval)
end = time.clock()


results = unique(nondominated(algorithm.result))



########### writing results

current_time = time.strftime("%H:%M:%S")
time_title = current_time[0:2]+'_'+current_time[3:5]+'_'+current_time[6:8]
date = time.strftime("%d/%m/%Y")
date_title = date[0:2]+'_'+date[3:5]+'_'+date[6:10]
computing_time = end-start

fileDir = os.path.dirname(os.path.realpath('__file__'))
filename = os.path.join(fileDir,'Results/'+date_title+'_'+time_title+'.txt')
 
file = open(filename,'w')
file.write('This optimization (simplified version 1) ended on '+date+' at '+current_time+'.\n') 
file.write('It lasted '+str(int(computing_time))+' seconds.\n')
file.write(str(class_for_opt.count)+' evaluations and '+str(nbr_year)+' years were used.\n') 
file.write('It provided '+str(len(algorithm.result))+' results')
file.write(', but only '+str(len(results))+' solution(s) are remaining after')
file.write(' deleting the twins and the dominated ones.\n')
file.write('Parameters were '+str(parameters)+'\n')
file.write('Results are '+str(results))

 
file.close() 
#print(results)
#print(get_value(b,'turbined_release',150))

#tracker.print_diff()
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:46:30 2017

@author: lraso
"""

#%% importing required modules 
import numpy as np #use log, matrix, matrix operations, (multivariate) normal ditribution
import pandas as pd
#import some_optimization_package as short_name #variable and constraints

from scipy.interpolate import interp1d

# Small change to test github

# %% class Element

class Element:
    """Base class 
    
    Element is the base class on which other elements are created by inheritance
    A new element created from Element must have the following methods:
    * __init__
    * _init_sim_generic
    * _sim
    
    
    """
    

    def __init__(self,name,system=self,variables=[],**kwargs): 
        
        """Create _Element object
        
        
        Attributes
        ----------
        :name: Element name
        :system: System or Element
                The system or element to which the element or sub-element belongs
                
        :variables: list of variables 
                Each variable is a dictionary having the following keys
                    * name
                    * min
                    * max
                    * val
        :delta_t: lenght of simuation timestep, in seconds [default: 1 day]
        :T: periodicity, in number of timesteps [default: 1 year, i.e. 365 steps] 
        
        * To be tested
        
        """ 
       
        # test uniqueness of name in system
        if any([name==element.name for element in system.internal_elements]):
            raise Exception('Element name non unique') # use error() instead?
            
        self.name=name 
        self.system=system        
        self.variables=variables
        
        
        self.internal_elements=internal_elements if internal_elements in kwargs else []
        
        self.input=[]
        for input_ in inputs:        
            add_input(input_) 
        #self.elements_down=elements_down if elements_down in kwargs else []
        
        
        
        if system==self        
            self.delta_t=delta_t if delta_t in kwargs else 24*60*60 # seconds per day, i.e. daily time steps    
            self.T=T if T in kwargs # system period [in number of steps]
        else:
            self.delta_t=system.delta_t
                    
        
        #create a list of variables of all components in system  
        
        for variable in variables: system.variables.append(self.name+ '_' + variable['name']):
        # Add dictionary to find element and variable from name_element_name_variable    
        

        
        self.__repr__()
        print(' has been created')
 

    
    def __repr__(self):
        print('Element ' + self.__class__.__name__)
        print('Name: ' + self.name)
        print(var['name']) for var in variables


    
    def add_input(self,element,var_name=[]):   #TODO add variable directly 
        
        """Add input to the element
        
        Attributes
        ----------
        :element: input element 
        :var_name: name of input variable 
        
        * To be tested
        
        """
        
        #if only one variable, take first var
        if len(element.variables)==1:
            var_name=element.variables[0]['name'] #take first variable's name
        
        
        #this can be worked out bettr and requries some checks              
        for el in self.system.internal_elements: #{'element' : element.name,'variable': var_name }
            if el.name==element.name:
                for var in  element.variables
                    if var['name']=var_name
                        self.input=element.var
                    else:
                        error('variable not found')



    def find_var(self,name):
		"""find a variable (i.e. a dictionary) in list of variables by name
           * To be tested
	
		"""
		# from https://stackoverflow.com/questions/8653516/python-list-of-dictionaries-search	
		return ([var for var in self.variables if var['name'] == name][0])



    def find_var_t(self,name,t)
		"""find a variable value by name and t
           * To be tested
		"""
		# from https://stackoverflow.com/questions/8653516/python-list-of-dictionaries-search	
		var=find_var(name)
		val=var.['val'][t]
		return(val)



    def set_constraint(var_name,value,**kwargs)    
        """
        
        Attributes
        ---------
        var_name
        
        value: float
        
        Optional attributes
        ------------------
        :function: function
        
        :argument: function's argument variable name
        
        
        """
        
        # fixed value constraint
        #assert # or other?
        if any var_name==variable['name'] for variable in self.variables:
            self.variables[var_name]=value
            
        else:
            error('Variable '+ var_name + ' not found' )      
            
        
        

    def set_delta_t(self,delta_t):
        """Assign delta_t to all system elements (or to all subcomponents)
        
        * To be tested
        
        """
        
        for element in self.internal_elements:
            if element.internal_elements == []:
                element.delta_t = self.delta_t
            else
                element.delta_t=set_delta_t(delta_t)# or other elements that has a defined deltat
                #self.T=_Element.T #Deltat[:][i]=_Element.Deltat #extract deltat


     
    
    def create_input(self):
        """ Create a df of all inputs with the longest period and uniform delta_t 
        
        * Not to be used in OptVsSim
        
        """
        
        df=[]
        for element in self.internal_elements:
            if element.type isinstance(pd.Series):
                df=pd.concat([df,element.timeseries])#add timeseries to df ??? TODO: how?

        # TODO FROM PANDAS df_cut = df.select longest period with data for all timeseries
   
        # TODO FROM PANDAS df_interp=interpole(self.delta_t) # uniform delta_t
   
        return(df_input)
        
    
 
 
    def order_elements(self):
        """reorder list of internal elements to match correct order of simulation
        
        TODO NOT WORKING YET. ELEMENTS MUST BE ADDED IN THE RIGHT ORDER
         
        * Not to be used in OptVsSim
        
        """
        
        order_list = [None] * len(self.internal_elements)
              
#        if False: 
#        #TODO
#            while any value in order_list is not None:
#                for element in self.internal_elements:
#                    order_values = order_list[]
#                    if element.inputs == []:
#                        index_element=system.internal_element.index(element)
#                        order_value [index_element]  = 0
#                    elif 
#                    all element.inputs is not None:
#                        order_value = max (order_value[element.inputs.index]) + 1
                                        

    def set_initial_conditions(self,var_name,initial_value):
        """Set initial condition for var to initial_value
        
         * Not to be used in OptVsSim
        """
        var=self.find_var(var_name)
        var['init']=initial_value #TODO this must be corrected


    def init_sim(self,H):
        """initialize simulation for all elements
        
        Attributes
        ----------
        :H: length of the simulation, in time steps
        
        """
        
            #create empty lists for all variables
        for element in self.internal_elements:
            element._init_sim(H)
	
	
    def _init_sim(self,H):  
        """initialize element simulation
        
        """
        for var in element.variables:
            var['val'].append([None] * (H-1)) #create empty list
            ##TODO ADD check : var['init'] exists and it is a valid value
            #var['val'][0]= var['init'] #set initial conditions to var[init]


    def _sim(self,t):
        """one step simulation
        
        run all internal elements
        
        """
        for element in self.internal_elements:
            element._sim(t)
            
        
        
    def sim_H(self):
        """
        
        TODO: get H from data
        """        
                
        #self.order_elements()
        
        #df_input=create_input()
        #H = len (df_input)
        H=1000
        
        self.init_sim()
            
        for t in range(H-1):
            self._sim(t)
            
        # create dataframe of all variables    

   
   
   def GP_sim(self,parameters):
       
       #TODO 
       for element in internal_elements:
           if isinstance(element,ReleaseRules):
                   element.parameters = parameters #lambda x : radial_function(x,parameters)
               
       self.sim_H()
       
       #calculated indicateors J_i
       # J_i=something
       
       return(J_i)
        
        
 
   ## Add _Element to set up SDDP problem in Problem                            #B I did not consider this part
    def add_to_SDDP_problem(self,system):                                       #B I did not consider this part
                
        
        if not isinstance(self,acDischargeDecision):                            #B I did not consider this part
            for t in range(Problem.T):                                          #B I did not consider this part
                for Var in self.Vars:                                           #B I did not consider this part
                    Problem.SDDP_Problems[t].A[Var.Pos][Var.Pos]=1              #B I did not consider this part


System=Element #just a name change

      
        
class Reservoir(Element):
    """ Reservoir element
    
    Attributes
    ----------    
    :system: System
            The system to which the Reservoir belongs to
        
    :max_volume: float
            Reservoir maximum volume [:math:`m^3`]
    
    :releases: var 
            variable with constraints(?) (having min=[0] and max value)
    
    
    Additional attributes (Optional)
    -------------------------------
    :release_policy:
    
    :specific_evaporation: timeseries
    
    :losses:
    
    :min_volume:
    
    :inflow: _Element
            Infow to the Reservoir [:math:`m^3/s`]
            
    :release_constraints:
    
    :surface: float
            Reservoir surface [:math:`m^2/s`]    
    :S_h_curve:
    
    :V_h_curve:
    
    :releases_h_down_curve:
    
    """
    
    
    def __init__(self,name,system,max_volume=+inf,release_name='release',max_release=+inf,**kwargs):
        
        min_volume=min_volume if 'min_volume' in kwargs else 0              
        volume={'name':'volume','max':max_volume,'min'=min_volume}
        
        for release in releases
            release={'name':release,'min':0,'max'=max_release} 
        #TODO        HOW to handle releases
        
        Element.__init__(self,name,system,[volume,releases])
        
               

        # TODO: add only if existing
        self.inflow=self.add_input(inflow) if inflow in kwargs     
        
        
        #Evaporation
        specific_evaporation= 0 if specific_evaporation not in kwargs
        
        if type (specific_evaporation) == float:
            specific_evaporation= pd.Series('1/1/1990',specific_evaporation)
            
        self.specific_evaporation=repeat(specific_evaporation,1/1/1900,now)
        
        self.losses=losses  if losses in kwargs                                     #B same here
        
        #Curves 
        S_v=S_h(inverse(v_h)) #PSEUDOCODE           #B don't really get what are the variables here.
        
        self.S_v=S_v #S_h_curve=( V_h_curve ^-1(dot))
    
        if  self.specific_evaporation != 0  and self.S_v=None:
            warning('Reservoir surface required')
            
    
        self.release_rules = []
          
    
    def set_evaporation(self,evap_ts):
        """ create array with H values of evaporation
         made of repetition of T values, taken from the timeseries evap_ts
        """
        evap_elab=1
        self.evap_specific=evap_elab        
        return()   
       
        
        
        
    def __repr__(self):
        _Element.__repr__(self)
        print('Operational Volume: %s m^3' % self.volume)    #B attribute volume has never been assigned 
        
          
            
            
    def set_initial_volume(self,value):
        """Set initial condition by initial volume
        
        
        """
        Element.set_initial_conditions('volume',value)
         
                


    def _sim(self,t):
        """
        
        """
        #TODO Add constraints on releases
        #TODO add constraint on volume
               
        x_t_min_1 = self.var['volume'][t]
        self.release_rules._sim(t)
        evaporation = self.find_var_t('evap',t)*self.S_v(x_t_min_1)
        #evap_t=self.find_var_t(evaporation) #is this correct?
        
        inflow = sum(self.inflow['val'][t]) #pseudocode  
        outflow = sum( release , evaporation , self.losses )
        
        self.var['release'][t] = releases
        self.var['volume'][t+1] = x_t_min_1 + self.delta_t * ( inflow - outflow )
        
        
		
		
	

    def add_to_SDDP_problem(self,Problem):
        #_Element.Add2Optim(Problem)        
        _Element.add_to_SDDP_problem(self,Problem)
        #presently, V_min and V_max cannot be timevariant
        
        for t in range(Problem.T):
            Problem.SDDP_Problems[t].LB[self.Vars[0].Pos]=self.minimum_volume #[t] if timevariant#B I did not consider this part
            Problem.SDDP_Problems[t].UB[self.Vars[0].Pos]=self.minimum_volume + self.volume      #B I did not consider this part
            Problem.SDDP_Problems[t].B[self.Vars[0].Pos][self.Vars[0].Pos]=1#no evaporation,by now#B I did not consider this part


## Inequality constraint. variable <= funct (arg) 

class ReleaseRules(Element):
    
    def __init__(self,name,reservoir,inputs,**kwargs):
        Element.__init__(self,name,system=reservoir,release)
    
    

# CLASS TO BE DELETED
class Constraint():
    def __init__(self,var,funct,arg,system,name=[]):

        self.smaller=smaller  #smaller is not defined yet
        self.name=name     
        
        self.add_constraint(self,system)
        
    
    
    def add_constraint(self,system);    
        
        #check variable presence         
#        for var in system.variables:
#            if variable is not in reservoir.variables:
#                error('variable ' + variable + ' not found' )
#                
        
        system.constraints.append(constraint)
#OR
    
 


## adding up discharge
class Confluence(_Element):
    """ confluence of two or more discharges
    
    """
    
    def __init__(self,name,system,elements_up):
         _Element.__init__(self,system,variables='discharge','elements_up'=elements_up)
         
    def simulate(self): 
        q_out=sum(something_here)
        return(q_out)



class Floodplain(_Element):
    """ Floodplain component
    Floodplain transforms the water-level time series  in flood surface
    """
    
    def __init__(self,name,system,S_h_table):
        _Element.__init__(self,system,{'name':'surface'})
        self.S_h_table=S_h_table        
        
    def simulate(self,h_ts):
        """
        
        
        """
                 
        ## TODO: see Bader's paper 
       
        return(surface_year)


class Inflow(Element,pd.Series):
    """Inflow is a timeseries.
    
    
    """
        
    def __init__(self,name,system,timeseries):
        
        variable={'name':'discharge'}
        _Element.__init__(self,system,variable)
        self.timeseries=timeseries
        self.order=0
    
    def _sim(self,t):
        pass

    def set_initial_conditions(self):
        pass

	def _init_sim(self,H):  #TODO: name to be modified!
		"""initialize simulation for timeseries
        
        Attributes
        ----------
        :H: length of the simulation [in time steps]
        
        * To be tested
        
        """                  
        var['val'] = self.timeseries.values[:H] #[0:H] it must be H values


class Indicator():
    """ Indicator of system performances
    
    Attributes
    ---------
    :input_var_name: str type
                    name input variable(s)
    
    :function: function
            function transforms the inout variable in the variable of interest
            example:  n of days for q> q^*, etc, is TODO:
            default value: do nothing /i.e. take input variable as it is
    
    :statistic: function
            function aggregating the variable of interestet
            example: 9th quantile 
            default: mean
    """
    
    def __init__(self,name,system,input_var_name,function=lambda x: x,statistic=mean):
       
        variable_in={'name':var_name,'I/O':'input'} # TODO make it generic for a list of n input_vars_name 
        variable_out={'name':name,'I/O':'input'}
        _Element.__init__(self,name,system,[variables_in,variable_out])

       #test input_var
        
        #TODO test if function has same number of inputs than n of input_vars_name
        self.function=function # function defining the indicator (examples: n of days for q> q^*, etc)
        self.statistic=statistic #statistic aggregating data
        
            #output   
        
    def evaluate(self):
        """ evaluate indicator
        
        """
        
        # TODO select inut and output type
        indicator_value=self.statistic(self.function(self.variables[input_var_name]))     
        return (indicator_value)
     
 
 
    
    def add_to_SDDP_problem(self,Problem): 
       #FIND pos=Problem._Elements.pos | self.Objectives(i)._Element_id==self._Elements.id
        VarPos= (Var.Pos for Var in self._Element.Vars if Var.Name==self.VarName)                
        for t in range(Problem.T):
            Problem.SDDP_Problems[t].f[VarPos]=self.W















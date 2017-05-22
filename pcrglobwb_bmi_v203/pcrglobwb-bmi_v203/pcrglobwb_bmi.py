# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 11:28:38 2014

@author: haag

Code heavily based on PCR-GLOBWB's own deterministic_runner.py and wflow's bmi 
(https://code.google.com/p/wflow/)
"""

# TO-DO: add function to get current timestep (not actual time)
# TO-DO: finish code for def finalize()
# TO-DO: make 'get_var' and 'set_var' work for all modules (currently "hardcoded" for:
#        groundwater, landSurface, landCover, meteo, pcrglobwb, routing, WaterBodies)

import logging
import numpy

from pcraster import *
from pcraster.framework import DynamicModel
from pcraster.framework import DynamicFramework

from deterministic_runner import DeterministicRunner
from configuration import Configuration
from currTimeStep import ModelTime
from spinUp import SpinUp

logger = logging.getLogger(__name__)

class pcrglobwbBMI(object):
    
    def initialize(self, config_file_location=None):
        """
        Initializes the model: read config file, load variables, get timestep information, etc.
        """
        self.configuration = Configuration(config_file_location)
        
        self.initial_state = None
		
        self.currTimeStep = ModelTime() # timeStep info: year, month, day, doy, hour, etc
        
        self.currTimeStep.getStartEndTimeSteps(self.configuration.globalOptions['startTime'],
                                               self.configuration.globalOptions['endTime'])
        
        self.deterministic_runner = DeterministicRunner(self.configuration, self.currTimeStep, self.initial_state)
        
        self.dynamic_framework = DynamicFramework(self.deterministic_runner,1)
        self.dynamic_framework.setQuiet(True)
        self.dynamic_framework._runInitial()
        self.dynamic_framework._runResume()
        
        # set timestep (to allow updates on a per-timestep-basis)
        self.currenttimestep = 0
		
        logger.info('Model initialized. Spin-up might be required.')
        
    def finalize(self):
        """
        Finalizes the model: shut down the model run, clean up resources, etc.
        """        
        self.dynamic_framework._runSuspend()
        #dynamic_framework._wf_shutdown()   # commented out, special function from wflow Dynamic Framework
		
    def spinup(self):
        """
        Spin-up the model. This is required to obtain realistic starting conditions for the model run.
        It runs on a yearly basis until the required convergence or max. allowed spin-up runs is reached.
        """        
        spin_up = SpinUp(self.configuration)                   # object for spin_up
        
        self.currTimeStep = ModelTime() # timeStep info: year, month, day, doy, hour, etc
        
        # spin-up
        noSpinUps = int(self.configuration.globalOptions['maxSpinUpsInYears'])
        if noSpinUps > 0:
            
            logger.info('Spin-Up #Total Years: '+str(noSpinUps))
            
            spinUpRun = 0 ; has_converged = False
            while spinUpRun < noSpinUps and has_converged == False:
                spinUpRun += 1
                self.currTimeStep.getStartEndTimeStepsForSpinUp(
                        self.configuration.globalOptions['startTime'],
                        spinUpRun, noSpinUps)
                logger.info('Spin-Up Run No. '+str(spinUpRun))
                deterministic_runner = DeterministicRunner(self.configuration, self.currTimeStep, self.initial_state)
                
                all_state_begin = deterministic_runner.model.getAllState() 
                
                self.dynamic_framework = DynamicFramework(deterministic_runner, self.currTimeStep.nrOfTimeSteps)
                self.dynamic_framework.setQuiet(True)
                self.dynamic_framework.run()
                
                all_state_end = deterministic_runner.model.getAllState() 
                
                has_converged = spin_up.checkConvergence(all_state_begin, all_state_end, spinUpRun, deterministic_runner.model.routing.cellArea)
                
                self.initial_state = deterministic_runner.model.getState()
		
		# setting model ready after spin-up
		self.currTimeStep.getStartEndTimeSteps(self.configuration.globalOptions['startTime'],
                                               self.configuration.globalOptions['endTime'])
		
		self.deterministic_runner = DeterministicRunner(self.configuration, self.currTimeStep, self.initial_state)
        
        logger.info('End of spin-up. Model is ready for transient simulation run.')
        
    def update(self, dt=-1):
        """
        Updates the model a number of timesteps, dependent on specified dt:
        dt = -1	-> runs the entire model from start time to end time
        dt = 1  -> updates the model 1 timestep (1 day)
        dt > 1  -> updates the model a number of timesteps (dt days)
        
        NOTE: the model can only run on a daily timestep!
        """
        
        if dt == 1:
            # update timestep
            self.currenttimestep += 1
            self.currTimeStep.update(self.currenttimestep)
            
            # commented out, already stated at initialization and at end of spin-up, not required at every timestep?
			#deterministic_runner = DeterministicRunner(self.configuration, self.currTimeStep, self.initial_state)
            
			# update model
            self.dynamic_framework = DynamicFramework(self.deterministic_runner, self.currenttimestep, self.currenttimestep)
            self.dynamic_framework.setQuiet(True)
            self.dynamic_framework.run()
            
            # update states (commented out, not required?)
            #self.initial_state = deterministic_runner.model.getState()

        elif dt == -1:
            # commented out, already stated at initialization and at end of spin-up, not required here as well?                                  
            #deterministic_runner = DeterministicRunner(self.configuration, self.currTimeStep, self.initial_state)
    
            self.dynamic_framework = DynamicFramework(self.deterministic_runner, self.currTimeStep.nrOfTimeSteps)
            self.dynamic_framework.setQuiet(True)
            self.dynamic_framework.run()
            
        else:
            # update timestep
            self.currenttimestep += 1
            self.currTimeStep.update(self.currenttimestep)
            
            # update model
            self.dynamic_framework = DynamicFramework(self.deterministic_runner, self.currenttimestep + (dt-1), self.currenttimestep)
            self.dynamic_framework.setQuiet(True)
            self.dynamic_framework.run()
            
            # update time
            self.currenttimestep += (dt-1)
            self.currTimeStep.update(self.currenttimestep)
    
    def get_start_time(self):
        """
        Returns model start time
        Input:  -
        Output: time as datetime (YYYY,MM,DD)
        """
        return self.currTimeStep.startTime
    
    def get_end_time(self):
        """
        Returns model end time
        Input:  -
        Output: time as datetime (YYYY,MM,DD)
        """
        return self.currTimeStep.endTime
    
    def get_current_time(self):
        """
        Returns current model time
        Input:  -
        Output: time as datetime (YYYY,MM,DD)
        """
        return self.currTimeStep.currTime
    
    def get_time_step(self):
        """
        Return current model timestep
        Input:  -
        Output: timestep as int
        """
        return self.currTimeStep.timeStepPCR

    def get_var(self, name, missingValues=-999):
        """
        Returns a numpy array from model library
        Input:  variable/map name (string)
        Output: numpy array or single variable, depending on input
        
        NOTE1: to get a variable from a specific landCover type, a tuple containing two strings should be used, with:
        - string 1 = name of landCover type
        - string 2 = name of variable
        
        NOTE2: there are two options to create a numpy array:
        - pcr2numpy    -> requires a value for MV (optional, default = -999)
        - pcr_as_numpy -> automatically sets nan for all MV
        Currently using pcr2numpy!
        """
        
        # check size of name input
        if numpy.size(name) == 1:
            
            # check for 'name' in the different sections of the model
            if hasattr(self.deterministic_runner.model.landSurface, name):
                pcrmap = getattr(self.deterministic_runner.model.landSurface, name)
            elif hasattr(self.deterministic_runner.model.routing, name):
                pcrmap = getattr(self.deterministic_runner.model.routing, name)
            elif hasattr(self.deterministic_runner.model.meteo, name):
                pcrmap = getattr(self.deterministic_runner.model.meteo, name)
            elif hasattr(self.deterministic_runner.model.groundwater, name):
                pcrmap = getattr(self.deterministic_runner.model.groundwater, name)
            else:
                logger.warn(name + " cannot be found in the model, returning empty list!")
        
        else:
            
            # first check if a specific model section was used as input
            if name[0] == 'landSurface':
                if hasattr(self.deterministic_runner.model.landSurface, name[1]):
                    pcrmap = getattr(self.deterministic_runner.model.landSurface, name[1])
            elif name[0] == 'routing':
                if hasattr(self.deterministic_runner.model.routing, name[1]):
                    pcrmap = getattr(self.deterministic_runner.model.routing, name[1])
            elif name[0] == 'WaterBodies':
                if hasattr(self.deterministic_runner.model.routing.WaterBodies, name[1]):
                    pcrmap = getattr(self.deterministic_runner.model.routing.WaterBodies, name[1])
            elif name[0] == 'pcrglobwb':
                if hasattr(self.deterministic_runner.model, name[1]):
                    pcrmap = getattr(self.deterministic_runner.model, name[1])
            # otherwise check if it is a variable from a landCover type
            else:
                # use the first entry of 'name' to find correct landCover type, second entry to find variable
                try:
                    if hasattr(self.deterministic_runner.model.landSurface.landCoverObj[name[0]], name[1]):
                        pcrmap = getattr(self.deterministic_runner.model.landSurface.landCoverObj[name[0]], name[1])
                    else:
                        logger.warn('(' + name[0] + ', ' + name[1] + ") cannot be found in the model, returning empty list!")
                except:
                    logger.warn('(' + name[0] + ', ' + name[1] + ") cannot be found in the model, returning empty list!")
            
        # attempt to create a numpy array, otherwise try to give the single value, or return empty list if this is both not possible
        try:
            return_value = pcr2numpy(pcrmap, missingValues)
            #return_value = pcr_as_numpy(pcrmap)
        except:
            try:
                return_value = pcrmap
            except:
                return []
        
        return return_value
        
    def set_var(self, name, var, missingValues=-999):
        """
        Sets a pcr map with values from a numpy array.
        Input:  variable/map name (string), values (numpy array or single value), missing values (optional, default = -999)
        Output: -
        
        NOTE: to set a variable from a specific landCover type, a tuple containing two strings should be used, with:
        - string 1 = name of landCover type
        - string 2 = name of variable
        """
        
        # try to create a pcr map from numpy array, otherwise just use the single value
        try:
            pcrmap = numpy2pcr(Scalar, var, missingValues)
        except:
            pcrmap = var
        
        # check if LDD (requires additional step)
        if 'lddMap' in name:
            pcrmap = ldd(pcrmap)
        
        # check size of name input
        if numpy.size(name) == 1:
            
            # update pcr map used in model with set values
            if hasattr(self.deterministic_runner.model.groundwater, name):
                exec "self.deterministic_runner.model.groundwater." + name + " = pcrmap"
            elif hasattr(self.deterministic_runner.model.landSurface, name):
                exec "self.deterministic_runner.model.landSurface." + name + " = pcrmap"
            elif hasattr(self.deterministic_runner.model.meteo, name):
                exec "self.deterministic_runner.model.meteo." + name + " = pcrmap"
            elif hasattr(self.deterministic_runner.model.routing, name):
                exec "self.deterministic_runner.model.routing." + name + " = pcrmap"
            else:
                logger.warn(name + " is not defined in the model, doing nothing!")
        
        else:
            
            # first check if a specific model section was used as input
            if name[0] == 'landSurface':
                if hasattr(self.deterministic_runner.model.landSurface, name[1]):
                    exec "self.deterministic_runner.model.landSurface." + name[1] + " = pcrmap"
            elif name[0] == 'routing':
                if hasattr(self.deterministic_runner.model.routing, name[1]):
                    exec "self.deterministic_runner.model.routing." + name[1] + " = pcrmap"
            elif name[0] == 'WaterBodies':
                if hasattr(self.deterministic_runner.model.routing.WaterBodies, name[1]):
                    exec "self.deterministic_runner.model.routing.WaterBodies." + name[1] + " = pcrmap"
            elif name[0] == 'pcrglobwb':
                if hasattr(self.deterministic_runner.model, name[1]):
                    exec "self.deterministic_runner.model." + name[1] + " = pcrmap"
            # otherwise check if it is a variable from a landCover type
            else:
                try:
                    if hasattr(self.deterministic_runner.model.landSurface.landCoverObj[name[0]], name[1]):
                        exec "self.deterministic_runner.model.landSurface.landCoverObj['" + name[0] + "']." + name[1] + " = pcrmap"
                    else:
                        logger.warn('(' + name[0] + ', ' + name[1] +  + ") is not defined in the model, doing nothing!")
                except:
                    logger.warn('(' + name[0] + ', ' + name[1] +  + ") is not defined in the model, doing nothing!")
    
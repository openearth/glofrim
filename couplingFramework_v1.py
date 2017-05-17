# -*- coding: utf-8 -*-
"""
Introduction:
-------------
This is the key script to couple PCR-GLOBWB with either Delft3D Flexible Mesh ("DFM") by Deltares or 
LISFLOOD-FP ("LFP") by University of Bristol.
The coupling is achieved by making use of the Basic Model Interface (BMI) which allows for initializing,
updating, data manipulation, and finalization of models from a shell-like environment as for instance this python-script. 
For couple Delft3D FM or LISFLOOD-FP, the python-module "bmi.wrapper" has to be loaded.
For Delft3D FM, any compiled version (>1.1.201.48898) has already implemented a BMI-compatible structure, and the 
required variables accessible to the user (downloadable upon request).
For LISFLOOD-FP, however, a specifically designed version needs to be compiled which is currently only available for
version 5.9 as the model is not generically BMI-compatible (also downloadable here).
Also for PCR-GLOBWB, a BMI-compatible version needs to be used (also downloadable here).

Running the script:
-------------------
To run the script, an set-file containing the required specifications and paths is necessary.
Using python, run this file along with the set-file as follows:
	python couplingFramework_v1.py default.set
	
Disclaimer:
-----------
The disclaimers of each component involved in this coupling (i.e. PCR-GLOBWB, LIFLOOD-FP, Delft3D Flexible Mesh, BMI Wrapper)
remain valid unless stated otherwise.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved in case of making use of this coupling script.
And of course don't forget to cite the associated paper as well :)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2017 Jannis Hoch

@author: Jannis Hoch M.Sc., Department of Physical Geography, Faculty of Geosciences, Utrecht University (j.m.hoch@uu.nl)
@date: 08-05-2017
"""

# -------------------------------------------------------------------------------------------------
# LOAD REQUIRED LIBRARIES
# -------------------------------------------------------------------------------------------------

import netCDF4
from distutils.util import strtobool
import pdb as pdb
import pylab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import pyproj as pyproj
import datetime
import bmi.wrapper
import pcrglobwb_203_30min_1way_prefactored as pcrglobwb_203_30min_1way_prefactored
from pcrglobwb_203_30min_1way_prefactored import pcrglobwb_bmi
from pcrglobwb_203_30min_1way_prefactored import disclaimer  
from coupling_PCR_FM import coupling_functions
from coupling_PCR_FM import model_functions
from coupling_PCR_FM import configuration
import time

# -------------------------------------------------------------------------------------------------
# IMPORT MODEL SETTINGS FROM INI-FILE
# -------------------------------------------------------------------------------------------------

config = configuration.Configuration()
config.parse_configuration_file(sys.argv[1])

# -------------------------------------------------------------------------------------------------
# SPECIFY MODEL SETTINGS
# -------------------------------------------------------------------------------------------------

model_type = config.model_type['model_type']                                                   

latlon = strtobool(config.general_settings['latlon'])
if latlon == False:
	inProj  = pyproj.Proj(init=config.model_settings['model_projection'])
use_Fluxes = strtobool(config.general_settings['use_Fluxes'])
use_RFS = strtobool(config.general_settings['use_RFS'])
verbose = strtobool(config.general_settings['verbose'])

# -------------------------------------------------------------------------------------------------
# SET PATHS TO .SO / .DLL FILES
# -------------------------------------------------------------------------------------------------

# these may be changed according to personal file and folder structure
if model_type == 'DFM':
    model_path 			= '/path/to/DFLOWFM/libdflowfm.so'
elif model_type == 'LFP':
    model_path 			= '/path/to/LISFLOODFP/liblisflood.so' 
else:
    sys.exit('\nno adequate model defined in ini-file - define either DFM or LFP!\n')

# -------------------------------------------------------------------------------------------------
# SPECIFY NUMERICAL SETTINGS
# -------------------------------------------------------------------------------------------------

nr_pcr_timesteps                      = int(config.numerical_settings['number_of_timesteps'])                      
update_step                           = int(config.numerical_settings['update_step'])  
                      
secPerDay                             = 86400.
end_time 							  = nr_pcr_timesteps * secPerDay
fraction_timestep 					  = secPerDay / update_step

threshold_inundated_depth             = float(config.numerical_settings['threshold_inundated_depth'])                         
threshold_inundated_depth_rivers      = float(config.numerical_settings['threshold_inundated_depth_rivers'])                         
threshold_inundated_depth_floodplains = float(config.numerical_settings['threshold_inundated_depth_floodplains'])                          

# other
missing_value_landmask                = 255                         # 05min resolution landmask
missing_value_pcr                     = -999

# -------------------------------------------------------------------------------------------------
# PLOT AND PRINT OPTIONS
# -------------------------------------------------------------------------------------------------

# set default figure size
pylab.rcParams['figure.figsize']      = (14.0, 7.0)
# set plot stuff for coloured plots of DFM water depths
my_cmap = matplotlib.cm.get_cmap('Blues_r')
my_cmap.set_under('seagreen')
my_cmap.set_bad('seagreen')

# -------------------------------------------------------------------------------------------------
# SET PATHS TO MODELS
# -------------------------------------------------------------------------------------------------

model_dir       	= config.model_settings['model_dir'] 
model_file      	= config.model_settings['model_file']
model_proj			= config.model_settings['model_projection']                                    

config_pcr       	=  config.PCR_settings['config_pcr']
landmask_pcr     	=  config.PCR_settings['landmask_pcr']
clone_pcr        	=  config.PCR_settings['clone_pcr']

# -------------------------------------------------------------------------------------------------
# INITIALIZE AND SPIN-UP PCR-GLOBWB
# -------------------------------------------------------------------------------------------------
                                  
# get start time of simulation
t_start = datetime.datetime.now()
# initiate logging and define folder for verbose-output
verbose_folder = model_functions.write2log(model_dir, model_file, latlon, use_Fluxes, use_RFS, t_start)

# print disclaimer
print '\n>>> Please consider reading the PCR-GLOBWB Disclaimer <<<\n' 
disclaimer.print_disclaimer()
time.sleep(5)

# initiate PCR-GLOBWB
model_pcr = pcrglobwb_203_30min_1way_prefactored.pcrglobwb_bmi.pcrglobwbBMI()
model_pcr.initialize(config_pcr)
print '\n>>> PCR-GLOBWB Initialized <<<\n' 

# spin-up PCR-GLOBWB
model_pcr.spinup()

# -------------------------------------------------------------------------------------------------
# INITIALIZING HYDRODYNAMIC MODEL
# -------------------------------------------------------------------------------------------------

# initiate hydraulic model
model_hydr = bmi.wrapper.BMIWrapper(engine = model_path, configfile = (os.path.join(model_dir, model_file)))
model_hydr.initialize()

print '\n>>> Hydrodynamic Model Initialized <<<\n' 

# -------------------------------------------------------------------------------------------------
# EXCTRACTING RELEVANT DATA FROM MODELS
# -------------------------------------------------------------------------------------------------

if model_type == 'DFM':
    
    #- retrieving data from Delft3D FM 
    x_coords, y_coords, z_coords, bottom_lvl, cell_points_fm, separator_1D, cellAreaSpherical, xz_coords, yz_coords, modelCoords, \
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FM(model_hydr, model_pcr, landmask_pcr, clone_pcr, use_RFS)
				
    print '\n>>> DFM data retrieved <<<\n'
         
elif model_type == 'LFP':
    
    #- retrieving data from LISFLOOD-FP
    dx, dy, DEM, bottom_lvl, H, waterDepth, rows, cols, \
                list_x_coords, list_y_coords, coupledFPindices, grid_dA, cellAreaSpherical, SGCQin, \
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FP(model_hydr, model_dir, model_pcr, landmask_pcr, clone_pcr, verbose_folder, use_RFS, verbose)

	# turning off separator between 1-D and 2-D which is only required for DFM but appears in subsequent function calls
    separator_1D = 0. 
    #- computing LFP-coordinates    
    modelCoords = coupling_functions.getVerticesFromMidPoints(list_x_coords, list_y_coords, dx, dy, verbose)
	
    print '\n>>> LFP data retrieved <<<\n'

#- computing PCR-coordinates
PCRcoords = coupling_functions.getPCRcoords(landmask_data_pcr)
		
# -------------------------------------------------------------------------------------------------
# COUPLING THE GRIDS
# -------------------------------------------------------------------------------------------------

# this is only required for plotting later, not for actual coupling process
CoupledCellsInfoAll = coupling_functions.coupleAllCells(modelCoords,PCRcoords)

# actual spatial coupling process, converting single indices of coupled PCR cells to double (array,column) indices
CoupleModel2PCR, CouplePCR2model, CoupledPCRcellIndices = coupling_functions.assignPCR2cells(landmask_pcr, modelCoords, verbose)

# saving plots of coupled cells to verbose-folder
if (verbose == True) and (use_RFS != True): 
    coupling_functions.plotGridfromCoords(PCRcoords, modelCoords)
    plt.savefig(os.path.join(verbose_folder , 'AllCells.png'))
    coupling_functions.plotGridfromCoords(CoupledCellsInfoAll[1],CoupledCellsInfoAll[0])
    plt.savefig(os.path.join(verbose_folder , 'CoupledCells.png'))   
    plt.close('all')

# -------------------------------------------------------------------------------------------------
# TURNING OFF CHANNELSTORAGE, WATERBODYSTORAGE, WATERBODIES AND RUNOFF TO CHANNELS
# -------------------------------------------------------------------------------------------------  

model_functions.noStorage(model_pcr, missing_value_pcr, CoupledPCRcellIndices, CouplePCR2model)

# -------------------------------------------------------------------------------------------------
# TURNING OFF ROUTING BY PCR IN COUPLED AREA
# -------------------------------------------------------------------------------------------------

model_functions.noLDD(model_pcr, CoupledPCRcellIndices, verbose_folder, verbose)

# -------------------------------------------------------------------------------------------------
# CALCULATE DELTA VOLUMES (DAY 1)
# first day outside loop, to make sure PCR time is at start time (timestep 1) and not at end of spin-up (timestep 365)
# ------------------------------------------------------------------------------------------------- 
 
# retrieving PCR-GLOBWB and converting it to m3/d
delta_volume_PCR_coupled = model_functions.calculateDeltaVolumes(model_pcr, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr)

# dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
delta_water_fm, verbose_volume = model_functions.calculateDeltaWater(CouplePCR2model, CoupleModel2PCR, delta_volume_PCR_coupled, cellAreaSpherical, fraction_timestep, model_type, use_Fluxes)

# saving PCR-GLOBWB output volumes and volumes used as input to hydraulic models to verbose-folder
if verbose == True:
	# initial file objects
    fo_PCR_V_tot = open(os.path.join(verbose_folder, 'delta_volume_PCR_coupled.txt'), 'w')
    fo_verbose_volume = open(os.path.join(verbose_folder, 'verbose_volume.txt'), 'w')
    # aggregate daily totals
    delta_volume_PCR_total_aggr = np.sum(delta_volume_PCR_coupled)
    verbose_volume_aggr = np.sum(verbose_volume)
    # write daily totals to file
    fo_PCR_V_tot.write(str(delta_volume_PCR_total_aggr) + os.linesep)
    fo_verbose_volume.write(str(verbose_volume_aggr) + os.linesep)

# check to ensure that volumes are equal, i.e. no errors in water balance
if np.round(np.sum(delta_volume_PCR_coupled),4) == np.round(np.sum(verbose_volume), 4):
    print 'PCR volume out: ', np.sum(delta_volume_PCR_coupled)
    print 'DFM volume in: ', np.sum(verbose_volume)
    pass
else:
    print 'PCR volume out: ', np.sum(delta_volume_PCR_coupled)
    print 'DFM volume in: ', np.sum(verbose_volume)
    sys.exit('\nDFM input volumes do not match PCR output volumes!\n')

# reshaping data for LISFLOOD-FP from list to arrays
if model_type == 'LFP':
    delta_water_fm = model_functions.fillFPgrid(model_hydr, coupledFPindices, delta_water_fm, DEM, verbose_folder, verbose)
  
# -------------------------------------------------------------------------------------------------
# FIRST UPDATE (DAY 1)
# -------------------------------------------------------------------------------------------------

print '\n>>> update 1 started <<<\n'

# updating arrays with computed additional volumes; array used depends on model specifications
if (model_type == 'LFP'):
    model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)

while model_hydr.get_current_time() < (model_pcr.get_time_step() * secPerDay):
    
    # updating DFM on user-specified time step
    if (model_type == 'DFM'):
        model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)
    
    # updating DFM or LFP on daily time step
    if (model_type == 'LFP'):
        model_hydr.update()

# ----------------------------------------------------------------------------------------------------
# UPDATE DFM FOR THE REST OF THE MODEL PERIOD
# ----------------------------------------------------------------------------------------------------

while model_pcr.get_time_step() < nr_pcr_timesteps:
    
    # retrieving PCR-GLOBWB and converting it to m3/d
    delta_volume_PCR_coupled = model_functions.calculateDeltaVolumes(model_pcr, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr)                                                                                                  
        
    # dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
    delta_water_fm, verbose_volume = model_functions.calculateDeltaWater(CouplePCR2model, CoupleModel2PCR, delta_volume_PCR_coupled, cellAreaSpherical, fraction_timestep, model_type, use_Fluxes, verbose)

    # saving PCR-GLOBWB output volumes and volumes used as input to hydraulic models to verbose-folder
    if verbose == True:
        # aggregate daily totals
        delta_volume_PCR_total_aggr = np.sum(delta_volume_PCR_coupled)
        verbose_volume_aggr = np.sum(verbose_volume)
        # write daily totals to file
        fo_PCR_V_tot.write(str(delta_volume_PCR_total_aggr) + os.linesep)
        fo_verbose_volume.write(str(verbose_volume_aggr) + os.linesep)
    
    # reshaping data for LISFLOOD-FP from list to arrays
    if model_type == 'LFP':
        delta_water_fm = model_functions.fillFPgrid(model_hydr, coupledFPindices, delta_water_fm, DEM, verbose_folder, verbose) 
    
    # updating arrays with computed additional volumes; array used depends on model specifications
    if (model_type == 'LFP'):
model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)      

    # update DFM unless it has has reached the same time as PCR
    while model_hydr.get_current_time() < (model_pcr.get_time_step() * secPerDay):
        
        # updating DFM on user-specified time step
        if (model_type == 'DFM'):
            model_functions.updateModel(model_hydr, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, model_type, verbose)      
        
        # updating DFM or LFP on daily time step
        if (model_type == 'LFP'):
            model_hydr.update()   
        
# ----------------------------------------------------------------------------------------------------
# END OF MODEL PERIOD REACHED
# ----------------------------------------------------------------------------------------------------
    
# get end time of simulation
t_end = datetime.datetime.now()
# update and finalize logging
model_functions.write2log(model_dir, model_file, latlon, use_Fluxes, use_RFS, t_start, t_end) 
# close files
if verbose == True:
    fo_PCR_V_tot.close()
    fo_verbose_volume.close()
	
model_hydr.finalize()

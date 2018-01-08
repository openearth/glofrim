# -*- coding: utf-8 -*-
"""
Attention:
-----------
THIS IS NOT YET WORKING!

Introduction:
-------------
This is the key script to couple PCR-GLOBWB with either Delft3D Flexible Mesh ("DFM") by Deltares or
LISFLOOD-FP ("LFP") by University of Bristol.
The coupling is achieved by making use of the Basic Model Interface (BMI) which allows for initializing,
updating, data manipulation, and finalization of models from a shell-like environment.
For couple Delft3D FM or LISFLOOD-FP, the python-module "bmi.wrapper" has to be loaded.
For Delft3D FM, any compiled version (>1.1.201.48898) has already implemented a BMI-compatible structure, and the
required variables accessible to the user.
For LISFLOOD-FP, however, a specifically designed version needs to be compiled which is currently only available for
version 5.9 as the model is not originally BMI-compatible.
Also for PCR-GLOBWB, a BMI-compatible version needs to be used which is provided in the downloaded zip-file.

Before running the framework, the two python code packages "pcrglobwb_bmi_v203" and "coupling_PCR_FM" need to be installed
by executing "python setup.py develop" in the respective folders.

Literature and sources:
-----------------------
	BMI         -> https://csdms.colorado.edu/wiki/BMI_Description
				-> http://www.sciencedirect.com/science/article/pii/S0098300412001252
	bmi.wrapper -> https://github.com/openearth/bmi-python

Running the script:
-------------------
To run the script, an ini-file containing the required specifications and paths is necessary.
Using python, run this file along with the ini/set-file and the env-file as follows:
	python couplingFramework_v1.1.py default.set paths.env

Disclaimer:
-----------
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

In addition to this disclaimer, the disclaimers of ALL models and code involved in this coupling remain valid.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved when using this coupling script.

Copyright (C) 2018 Jannis Hoch

@author: Jannis Hoch, Department of Physical Geography, Faculty of Geosciences, Utrecht University (j.m.hoch@uu.nl)
@date: 08-01-2018
"""

# -------------------------------------------------------------------------------------------------
# LOAD REQUIRED LIBRARIES
# -------------------------------------------------------------------------------------------------

import netCDF4
from distutils.util import strtobool
import matplotlib
import matplotlib.pyplot as plt
import sys
import os
import pdb
import platform
import numpy as np
import pyproj as pyproj
import datetime
import time
import bmi.wrapper
import pcrglobwb_bmi_v203 as pcrglobwb_bmi_v203
from pcrglobwb_bmi_v203 import pcrglobwb_bmi
from pcrglobwb_bmi_v203 import disclaimer
from coupling_PCR_FM import coupling_functions
from coupling_PCR_FM import model_functions
from coupling_PCR_FM import configuration
from coupling_PCR_FM import utils as utils

# -------------------------------------------------------------------------------------------------
# IMPORT MODEL SETTINGS FROM INI-FILE
# -------------------------------------------------------------------------------------------------

# parse set/ini-file with central/general settings for coupling framework
config = configuration.Configuration()
config.parse_configuration_file(sys.argv[1])

# parse env-file for user-specific paths and environmental variables
envs = configuration.Configuration()
envs.parse_configuration_file(sys.argv[2])

# -------------------------------------------------------------------------------------------------
# SPECIFY MODEL SETTINGS
# -------------------------------------------------------------------------------------------------

type_hydrodynamicModel = config.hydrodynamic_model['model_name']
type_routingModel          = config.routing_model['model_name']
type_hydrologicModel      = config.hydrologic_model['model_name']

latlon = strtobool(config.general_settings['latlon'])
if latlon == False:
	inProj  = pyproj.Proj(init=config.hydrodynamic_model['model_projection'])

use_Fluxes = strtobool(config.general_settings['use_Fluxes'])
use_RFS = strtobool(config.general_settings['use_RFS'])
verbose = strtobool(config.general_settings['verbose'])

# -------------------------------------------------------------------------------------------------
# SPECIFY NUMERICAL SETTINGS
# -------------------------------------------------------------------------------------------------

nr_timesteps                      = int(config.numerical_settings['number_of_timesteps'])
update_step                           = int(config.numerical_settings['update_step'])

secPerDay                             = 86400.
end_time 							  = nr_timesteps * secPerDay
fraction_timestep 					  = secPerDay / update_step

threshold_inundated_depth             = float(config.numerical_settings['threshold_inundated_depth'])
threshold_inundated_depth_rivers      = float(config.numerical_settings['threshold_inundated_depth_rivers'])
threshold_inundated_depth_floodplains = float(config.numerical_settings['threshold_inundated_depth_floodplains'])

# other
missing_value_landmask                = 255
missing_value_pcr                     = -999

# -------------------------------------------------------------------------------------------------
# GET ENVIRONMENT SETTINGS
# -------------------------------------------------------------------------------------------------

inputDIR_env 	= envs.PCRpaths['inputDirectoryPCR']
outputDIR_env 	= envs.PCRpaths['outputDirectoryPCR']
DFM_path		= envs.DFM_engine['DFM_path'] #absolute
CMF_path		= envs.CMF_engine['CMF_path'] #relative

# -------------------------------------------------------------------------------------------------
# SET PATHS TO MODELS
# -------------------------------------------------------------------------------------------------

# get main folder (i.e. GitHub folder)
mainFolder = os.getcwd()

# hydrodynamic model
hydrodynamicModel_dir       	= config.hydrodynamic_model['model_dir']
hydrodynamicModel_dir           = os.path.join(mainFolder, hydrodynamicModel_dir)
hydrodynamicModel_file      	= config.hydrodynamic_model['model_file']
hydrodynamicModel_proj		    = config.hydrodynamic_model['model_projection']

# routing model
routingModel_dir       			= config.routing_model['model_dir']
routingModel_dir           		= os.path.join(mainFolder, routingModel_dir)
routingModel_file      			= config.routing_model['model_file']

# hydrologic model
hydrologicModel_dir				= config.hydrologic_model['config_dir']
hydrologicModel_file       		= config.hydrologic_model['config_file']
hydrologicModel_file_tmp		= str('tmp_'+hydrologicModel_file)
hydrologicModel_config	        = os.path.join(mainFolder, os.path.join(hydrologicModel_dir, hydrologicModel_file))
hydrologicModel_config_tmp 		= os.path.join(mainFolder, os.path.join(hydrologicModel_dir, hydrologicModel_file_tmp))

# overwriting inputDIR and outputDIR in PCR-GLOBWB ini-file with respect to personal environment
kwargs = dict(inputDir = inputDIR_env, outputDir = outputDIR_env)
utils.write_ini(hydrologicModel_config_tmp, hydrologicModel_config, **kwargs)

# parsing the overwritten PCR-GLOBWB ini-file
configPCR 						= configuration.Configuration()
configPCR.parse_configuration_file(hydrologicModel_config_tmp)

# retrieving informatino directly from PCR-GLOBWB ini-file
inputDIR 						= configPCR.globalOptions['inputDir']
clone_pcr 						= os.path.join(inputDIR, configPCR.globalOptions['cloneMap'])
landmask_pcr 					= os.path.join(inputDIR, configPCR.globalOptions['landmask'])

# -------------------------------------------------------------------------------------------------
# SET PATHS TO .SO / .DLL FILES
# -------------------------------------------------------------------------------------------------

if type_hydrodynamicModel == 'DFM':
    if platform.system() == 'Linux':
        path_to_hydrodynamicModel = DFM_path  # for Linux
    elif platform.system() == 'Windows':
         path_to_hydrodynamicModel = '/path/to/DFLOWFM/lib/libdflowfm.dll'  # for Windows

elif type_hydrodynamicModel == 'LFP':
    if platform.system() == 'Linux':
        path_to_routingModel = './lisflood-bmi-v5.9/liblisflood.so'  # for Linux
    elif platform.system() == 'Windows':
        sys.exit('\nLFP v5.9 with BMI currently not supported on Windows!\n')

if type_routingModel == 'CMF':
	if platform.system() == 'Linux':
		path_to_routingModel = os.path.join(mainFolder, CMF_path)
	elif platform.system() == 'Windows':
		sys.exit('\nCaMa-Flood with BMI currently not supported on Windows!\n')
elif type_routingModel != 'CMF':
	sys.exit('\nNo supported routing model chosen in ini/set-file!')
#TODO: somewhere in the model functions (think where geometry is extracted) I have a similar exit-statement; once
# CMF functions are written, add such statement there instead of here


# -------------------------------------------------------------------------------------------------
# INITIALIZE AND SPIN-UP HYDROLOGIC MODEL
# -------------------------------------------------------------------------------------------------

# get start time of simulation
t_start = datetime.datetime.now()
# initiate logging and define folder for verbose-output
#TODO: with new models being added (esp. new type of routing model) this functions needs an update
verbose_folder = model_functions.write2log(hydrodynamicModel_dir, hydrodynamicModel_file, latlon, use_Fluxes, use_RFS, t_start)

# print disclaimer
print '\n>>> Please consider reading the PCR-GLOBWB Disclaimer <<<\n'
disclaimer.print_disclaimer()
time.sleep(2)

# initiate PCR-GLOBWB
hydrologicModel = pcrglobwb_bmi_v203.pcrglobwb_bmi.pcrglobwbBMI()
hydrologicModel.initialize(hydrologicModel_config_tmp)
print '\n>>> Hydrologic Model Initialized <<<\n'

# spin-up PCR-GLOBWB
hydrologicModel.spinup()

# -------------------------------------------------------------------------------------------------
# INITIALIZING ROUTING MODEL
# -------------------------------------------------------------------------------------------------

routingModel = bmi.wrapper.BMIWrapper(engine = path_to_routingModel, configfile = (os.path.join(routingModel_dir, routingModel_file)))
routingModel.initialize()
print '\n>>> Routing Model Initialized <<<\n'

# -------------------------------------------------------------------------------------------------
# INITIALIZING HYDRODYNAMIC MODEL
# -------------------------------------------------------------------------------------------------

hydrodynamicModel = bmi.wrapper.BMIWrapper(engine = path_to_hydrodynamicModel, configfile = (os.path.join(hydrodynamicModel_dir, hydrodynamicModel_file)))
hydrodynamicModel.initialize()
print '\n>>> Hydrodynamic Model Initialized <<<\n'

# -------------------------------------------------------------------------------------------------
# EXCTRACTING RELEVANT DATA FROM MODELS
# -------------------------------------------------------------------------------------------------

if type_hydrodynamicModel == 'DFM':

    #- retrieving data from Delft3D FM
    x_coords, y_coords, z_coords, bottom_lvl, cell_points_fm, separator_1D, cellAreaSpherical, xz_coords, yz_coords, modelCoords, \
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FM(hydrodynamicModel, hydrologicModel, landmask_pcr, clone_pcr, use_RFS)
    print '\n>>> DFM data retrieved <<<\n'

elif type_hydrodynamicModel == 'LFP':

    #- retrieving data from LISFLOOD-FP
    dx, dy, DEM, bottom_lvl, H, waterDepth, rows, cols, \
                list_x_coords, list_y_coords, coupledFPindices, grid_dA, cellAreaSpherical, SGCQin, \
                cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_FP(hydrodynamicModel, hydrodynamicModel_dir, hydrologicModel, landmask_pcr, clone_pcr, verbose_folder, use_RFS, verbose)

    separator_1D = 0. # setting separator between 1-D and 2-D to 0 as only used for DFM

    #- computing FP-coordinates
    modelCoords = coupling_functions.getVerticesFromMidPoints(list_x_coords, list_y_coords, dx, dy, verbose)
    print '\n>>> LFP data retrieved <<<\n'

#- computing PCR-coordinates
PCRcoords = coupling_functions.getPCRcoords(landmask_data_pcr)

# -------------------------------------------------------------------------------------------------
# COUPLING THE GRIDS
# -------------------------------------------------------------------------------------------------

# this is only required for plotting later, not for actual coupling process
CoupledCellsInfoAll = coupling_functions.coupleAllCells(modelCoords,PCRcoords)

# converting single indices of coupled PCR cells to double (array,column) indices
CoupleModel2PCR, CouplePCR2model, CoupledPCRcellIndices = coupling_functions.assignPCR2cells(landmask_pcr, modelCoords, verbose)

# saving plots of coupled cells to verbose-folder
# currently doesn't work with FM and use_RFS on, due to data structure required (? check this ?)
if verbose == True:
    coupling_functions.plotGridfromCoords(PCRcoords, modelCoords)
    plt.savefig(os.path.join(verbose_folder , 'AllCells.png'))
    coupling_functions.plotGridfromCoords(CoupledCellsInfoAll[1],CoupledCellsInfoAll[0])
    plt.savefig(os.path.join(verbose_folder , 'CoupledCells.png'))
    plt.close('all')

# -------------------------------------------------------------------------------------------------
# TURNING OFF CHANNELSTORAGE, WATERBODYSTORAGE, WATERBODIES AND RUNOFF TO CHANNELS
# -------------------------------------------------------------------------------------------------

model_functions.noStorage(hydrologicModel, missing_value_pcr, CoupledPCRcellIndices, CouplePCR2model)

# -------------------------------------------------------------------------------------------------
# TURNING OFF ROUTING BY PCR IN COUPLED AREA
# -------------------------------------------------------------------------------------------------

model_functions.noLDD(hydrologicModel, CoupledPCRcellIndices, verbose_folder, verbose)

# -------------------------------------------------------------------------------------------------
# CALCULATE DELTA VOLUMES (DAY 1)
# first day outside loop, to make sure PCR time is at start time (timestep 1) and not at end of spin-up (timestep 365)
# -------------------------------------------------------------------------------------------------

# retrieving PCR-GLOBWB and converting it to m3/d
delta_volume_PCR_coupled = model_functions.calculateDeltaVolumes(hydrologicModel, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr)

# dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
delta_water_fm, verbose_volume = model_functions.calculateDeltaWater(CouplePCR2model, CoupleModel2PCR, delta_volume_PCR_coupled, cellAreaSpherical, fraction_timestep, type_hydrodynamicModel, use_Fluxes)

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
if np.round(np.sum(np.asfarray(delta_volume_PCR_coupled)), 3) == np.round(np.sum(np.asfarray(verbose_volume)), 3):
    print 'PCR volume out: ', np.sum(np.asfarray(delta_volume_PCR_coupled))
    print 'FM volume in: ', np.sum(np.asfarray(verbose_volume))
    pass
else:
    print 'PCR volume out: ', np.sum(np.asfarray(delta_volume_PCR_coupled))
    print 'FM volume in: ', np.sum(np.asfarray(verbose_volume))
    sys.exit('\nFM input volumes do not match PCR output volumes!\n')

# reshaping data for LISFLOOD-FP from list to arrays
if type_hydrodynamicModel == 'LFP':
    delta_water_fm = model_functions.fillLFPgrid(hydrodynamicModel, coupledFPindices, delta_water_fm, DEM, verbose_folder, verbose)

# -------------------------------------------------------------------------------------------------
# FIRST UPDATE (DAY 1)
# -------------------------------------------------------------------------------------------------

print '\n>>> update 1 started <<<\n'

# updating arrays with computed additional volumes; array used depends on model specifications
if (type_hydrodynamicModel == 'LFP'):
    model_functions.updateModel(hydrodynamicModel, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, type_hydrodynamicModel, verbose)

while hydrodynamicModel.get_current_time() < (hydrologicModel.get_time_step() * secPerDay):

    # updating FM on user-specified time step
    if (type_hydrodynamicModel == 'DFM'):
        model_functions.updateModel(hydrodynamicModel, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, type_hydrodynamicModel, verbose)

    # updating FM or FP on daily time step
    if (type_hydrodynamicModel == 'LFP'):
        hydrodynamicModel.update()

# ----------------------------------------------------------------------------------------------------
# UPDATE FM FOR THE REST OF THE MODEL PERIOD
# ----------------------------------------------------------------------------------------------------

while hydrologicModel.get_time_step() < nr_timesteps:

    # retrieving PCR-GLOBWB and converting it to m3/d
    delta_volume_PCR_coupled = model_functions.calculateDeltaVolumes(hydrologicModel, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr)

    # dividing delta volume from PCR-GLOBWB over hydraulic cells, depending on model specifications
    delta_water_fm, verbose_volume = model_functions.calculateDeltaWater(CouplePCR2model, CoupleModel2PCR, delta_volume_PCR_coupled, cellAreaSpherical, fraction_timestep, type_hydrodynamicModel, use_Fluxes)

    # saving PCR-GLOBWB output volumes and volumes used as input to hydraulic models to verbose-folder
    if verbose == True:
        # aggregate daily totals
        delta_volume_PCR_total_aggr = np.sum(delta_volume_PCR_coupled)
        verbose_volume_aggr = np.sum(verbose_volume)
        # write daily totals to file
        fo_PCR_V_tot.write(str(delta_volume_PCR_total_aggr) + os.linesep)
        fo_verbose_volume.write(str(verbose_volume_aggr) + os.linesep)

    # reshaping data for LISFLOOD-FP from list to arrays
    if type_hydrodynamicModel == 'LFP':
        delta_water_fm = model_functions.fillLFPgrid(hydrodynamicModel, coupledFPindices, delta_water_fm, DEM, verbose_folder, verbose)

    # updating arrays with computed additional volumes; array used depends on model specifications
    if (type_hydrodynamicModel == 'LFP'):
        model_functions.updateModel(hydrodynamicModel, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, type_hydrodynamicModel, verbose)

    # update FM unless it has has reached the same time as PCR
    while hydrodynamicModel.get_current_time() < (hydrologicModel.get_time_step() * secPerDay):

        # updating FM on user-specified time step
        if (type_hydrodynamicModel == 'DFM'):
            model_functions.updateModel(hydrodynamicModel, delta_water_fm, update_step, separator_1D, use_Fluxes, use_RFS, type_hydrodynamicModel, verbose)

        # updating FM or FP on daily time step
        if (type_hydrodynamicModel == 'LFP'):
            hydrodynamicModel.update()

# ----------------------------------------------------------------------------------------------------
# END OF MODEL PERIOD REACHED
# ----------------------------------------------------------------------------------------------------

# get end time of simulation
t_end = datetime.datetime.now()
# update and finalize logging
model_functions.write2log(hydrodynamicModel_dir, hydrodynamicModel_file, latlon, use_Fluxes, use_RFS, t_start, t_end)
# close files
if verbose == True:
    fo_PCR_V_tot.close()
    fo_verbose_volume.close()

#- finalizing hydrodynamic model to properly end execution
hydrodynamicModel.finalize()

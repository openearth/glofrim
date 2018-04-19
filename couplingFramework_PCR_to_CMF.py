
# ---
# LOAD lIBS
# ---

# load general libraries
import netCDF4
from netCDF4 import Dataset
import rasterio
import os, sys
from datetime import datetime
import numpy as np
import spotpy
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
# load local libraries
from coupling_PCR_FM.model_functions_v2 import PCR_model, CMF_model, DFM_model
from coupling_PCR_FM import coupling_functions_v2
from coupling_PCR_FM.utils import config_to_dict, determineSteps

# ---
# IMPORT MODEL SETTINGS FROM INI-FILE
# ---

argv1 = sys.argv[1] # settings file
argv2 = sys.argv[2] # paths file

# ---
# PARSE CONTENT
# ---

# parse set/ini-file with central/general settings for coupling framework
config = config_to_dict(argv1)
# parse env-file for user-specific paths and environmental variables
envs = config_to_dict(argv2)
# combine
config.update(envs)
options = config
# parse dates
start_date = datetime.strptime(options['numerical_settings']['startTime'], '%Y-%m-%d')
end_date = datetime.strptime(options['numerical_settings']['endTime'], '%Y-%m-%d')
timeSteps = determineSteps(start_date, end_date)

# ---
# SET OUTPUT DIRs
# ---

cwd = os.getcwd() # note: this get changed by pcr initialization later on
out_dir = options['PCRpaths']['outputDirectoryPCR']
out_dir = out_dir + 'PCR2CMF2DFM/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# ---
# CREATE PCR BMI OBJECT
# ---

PCR_config_fn = os.path.join(cwd, options['hydrologic_model']['config_dir'], options['hydrologic_model']['config_file'])
PCR_in_dir = options['PCRpaths']['inputDirectoryPCR']
PCR_out_dir = os.path.join(out_dir, 'PCR')

PCR_bmi = PCR_model(PCR_config_fn, PCR_in_dir, PCR_out_dir,
                            start_date, end_date)

# ---
# CREATE CMF BMI OBJECT
# ---

CMF_engine = os.path.join(cwd, options['CMF_engine']['CMF_path'])
CMF_model_dir = os.path.join(cwd, options['routing_model']['model_dir'])
CMF_config_fn = os.path.join(CMF_model_dir, options['routing_model']['model_file'])
CMF_out_dir = os.path.join(out_dir, 'CMF')

CMF_bmi = CMF_model(CMF_engine, CMF_config_fn, CMF_model_dir, CMF_out_dir,
                         start_date, end_date, dt=86400)

# ---
# COUPLE PCR TO CMF
# ---

PCR_bmi.couple_grid_to_grid(CMF_bmi)

# ---
# INITIALIZE CMF
# ---

CMF_bmi.initialize()

# ---
# INITIALIZE PCR
# ---

PCR_bmi.initialize()

# ---
# RUN COUPLED PCR->CMF MODEL
# ---

tStart = datetime.now()
for i in range(timeSteps):
    PCR_bmi.update()
    coupling_functions_v2.set_CMF_forcing(PCR_bmi, CMF_bmi)
    CMF_bmi.update()
tEnd = datetime.now()

print 'start time coupling: ', tStart
print 'end time coupling: ', tEnd
print 'average time per update PCR->CMF->DFM: ', abs((tEnd - tStart)) / timeSteps

# ---
# FINALIZE MODELS
# ---

PCR_bmi.finalize()
CMF_bmi.finalize()

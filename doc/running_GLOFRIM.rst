.. _running_GLOFRIM:

***************
Running GLOFRIM
***************

GLOFRIM exists of a series of uniformed BMI wrappers for each model and a overarching BMI wrapper for running coupled models.

Run GLOFRIM from Python
=======================
Coupled run
-----------

To run a coupled model from python use the following lines. 
The glofrim.ini (see example in root directory) configuration file hold the information of the individual model configuration files and exchanges between the models::

  # import GLOFRIM
  from glofrim import Glofrim 
  # initialize coupled bmi
  cbmi = Glofrim() 
  # initialize the coupling with the glofrim.ini configuration file
  cbmi.initialize_config(/path/to/glofrim.ini) 

A basic model run uses the following statements::

  # optional: get the model start time
  bmi.get_start_time() 
  # initialize model
  bmi.initialize_model() 
  # run until set endtime
  bmi.update_until(bmi.get_end_time()) 
  # finalize model
  bmi.finalize()

Stand-alone run
---------------

To run stand alone models via the GLOFRIM BMI wrapper you can use the lines below, followed by the same statements as before::

  # import the CaMa-Flood bmi wrapper
  from glofrim import CMF 
  # intialize bmi with reference to engine (only for non-python models)
  bmi = CMF(/path/to/model_engine) 
  bmi.initialize_config(/path/to/model_configuration_file)


Run GLOFRIM from command line
=============================

The GLOFRIM library contains a script to run combined and single models with a single line from a terminal. 
This script can be found in the glofirm-py/scripts folder.

GLOFRIM can be executed as follows on Linux command line::

  python glofrim_runner.py run /path/to/glofrim.ini --env /path/to/glofrim.env -s startdate -e enddate

Both *startdate* and *enddate* must be in yyyy-mm-dd format.

For more info on coupled runs, check::

  python glofrim_runner.py run –help

and for stand-alone runs::

  python glofrim_runner.py run_single –help

.. _the_ini_file:

The GLOFRIM configuration file
==============================
The GLOFRIM configuration (or .ini) file has four sections, the engines, models, coupling and exchanges settings, each is shortly explained here.


The **engines** section contains paths to the shared libraries of each non-python model. For convenience the absolute paths in the engine and models sections 
may also be set in a seperate environment.env file in the GLOFRIM root folder.::

    [engines]
    # path to model engines; only required for the non-python models used
    # these settings can also be set in environment.env
    CMF = /path/to/libcama.so
    DFM = /path/to/libdflowfm.so
    LFP = /path/to/liblisflood.so

The **models** section needs the paths to all model configuraiton files. Together with the model engine, this allows GLOFRIM to know the model schematisation and to 
communicate with that model via BMI. Add only models which are part of the (coupled) run. The paths should be either relative to the root_dir option (if set), this ini file directory or absolute.::

    [models]
    # alternative root dir for relative ini-file paths, by default the directory of this ini file is used; 
    # this setting can also be set in environment.env
    root_dir = /path/to/models

    # all models which are listed here are run during update
    # format: model_short_name = /path/to/configuration_file 
    PCR=/path/to/pcrglobwb.ini
    WFL=/path/to/wflow.ini
    CMF=../rel_path/to/input_flood.nam.org
    LFP=/path/to/lisflood.par
    DFM=rel_path/to/dflowfm.mdu

.. note::
    Note that the user can change model options through the GLOFRIM API. For all models but WFL, a new configuration file name ending with *_glofrim* is written to communicate these changes with the model
    before model initialization. For WFL it's possible to communicate these changes directly via BMI. 

.. note::
    CMF only listens to the configuration file if it is called *input_flood.nam*, therefore the original configuration file should be called different, for instance input_flood.nam.org.

The **coupling** section contains general settings for the exchanges between models.::

    [coupling]
    # timestep for exchanges [sec]
    dt=86400

The **exchanges** section contains the information about how the models communicate on run time. This part has a slightly complex syntax as it contains a lot of information. 
Every line indicates one exchange from the left (upstream/get) model.variable to the right (downstream/set) model.variable. This can be further extended by multipliers which can be model variables 
or scalar values in order to make sure the variable units match. Finally behind the @ the spatial location to get (upstream) and set (downstream) the model variables. 
Current options are @1d,  @1d_us (the most upstream 1d cells or nodes) and @grid_us (the upstream cell for each grid cell)::

    [exchanges]
    # setup exchanges which are executed during the coupled update function. 
    # format: From_model.var1*var2*multiplier@index = To_model.var*multiplier@index
    # the multiplier is optional; if no index is set, by default the whole 2D domain is coupled

    # Example 1: PCR runoff [m] to CMF runoff [m] 
    # The interal CMF interpolation matrix is used to convert from the PCR grid to the CMF U-Grid.
    PCR.runoff=CMF.roffin 

    # Example 2: PCR runoff [m] & upstream discharge [m3/s] to DFM rain [mm] (used as api for lateral inflows) 
    # both sides are converted to volumes per exchange timestep [m3/day]
    PCR.runoff*cellArea=DFM.rain*ba*1000@1d
    PCR.discharge*86400@grid_us=DFM.rain*ba*1000@1d_us

.. note::
    Note that only exchange of fluxes has been tested so far. 


.. _running_GLOFRIM:


***************
Running GLOFRIM
***************

GLOFRIM exists of a series of uniformed BMI wrappers for each model and a overarching BMI wrapper for running coupled models.

Coupled run
===========

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
===============

To run stand alone models via the GLOFRIM BMI wrapper you can use followed by the same statements as before::

  # import the CaMa-Flood bmi wrapper
  from glofrim import CMF 
  # intialize bmi with reference to engine (only for non-python models)
  bmi = CMF(/path/to/model_engine) 
  bmi.initialize_config(/path/to/model_configuration_file)

Convenience script
==================
The GLOFRIM library contains a script to run combined and single (for testing purposes) models with a single line from a terminal. 
This script can be found in the glofirm-py/scripts folder.

GLOFRIM can be executed as follows on Linux command line::

  python glofrim_runner.py run /path/to/glofrim.ini --env /path/to/glofrim.env -s startdate -e enddate

Both *startdate* and *enddate* must be in yyyy-mm-dd format.

Help
====

For more info on coupled runs, check::

  python glofrim_runner.py run –help

and for stand-alone runs::

  python glofrim_runner.py run_single –help


   

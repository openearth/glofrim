# GLOFRIM 2.0
Globally Applicable Framework for Integrated Hydrological-Hydrodynamic Modelling (GLOFRIM)

development by Jannis M. Hoch (Utrecht University, Deltares), Dirk Eilander (VU Amsterdam, Deltares), and Hiroaki Ikeuchi (University of Tokyo) \
contact: Jannis M. Hoch (j.m.hoch@uu.nl), Dirk Eilander (dirk.eilander@vu.nl)

We also want to acknowledge the contributions of all colleagues at the insitues involved in the development of GLOFRIM.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597107.svg)](https://doi.org/10.5281/zenodo.597107)

# Description
GLOFRIM offers a flexible and modular tool to couple hydrologic, routing, and hydodynamic models across scales. This enables integration of physical processes from different models. The coupling process is spatially explicit (i.e. on grid-to-grid basis) and model information is exchanged online (i.e. per time step).

GLOFRIM is designed as a “human interface” with additional and user friendly Python functions on top of the basic model interface (BMI), which makes it easy to setup and run coupled model simulations. For the model developer, only the BMI needs to be implemented in the model in a scripting language of choice, which makes it easy to develop and maintain. 

While version 1 allowed for coupling PCR-GLOBWB with either Delft3D Flexible Mesh or LISFLOOD-FP, version 2 has a more generic setup and has been extended with the hydrologic modelling suite wflow and the global routing model CaMa-Flood.

With the available models, different coupled hydrologic and hydrodynamic model runs can be done, for instance:
 - 2-step coupling: hydrology -> 1D routing or hydrology -> full 2D hydrodynamics
 - 3-step coupling: hydrology -> 1D routing -> full 2D hydrodynamics

![alt text](/doc/_images/GLOFRIM_flows_wLegend.png "Conceptual GLOFRIM diagram")

## Warrants
Currently, the coupling process in only one-directional, i.e. only downstream along the model cascade.
Work is currently performed to extend it to a full feedback loop.

It is important to note that GLOFRIM provides only a tool to coupled models across scales and processes. The quality of simulations therefore still depends on the quality of the model discretizations used.

The framework has successfully been tested on Linux platforms. 
Please note that it running it on Windows is currently not supported.

## Model specs
 - PCR-GLOBWB: since the model does not generically contain BMI function, a bespoke version is provided with the GLOFRIM package.
 - Delft3D Flexible Mesh: the model is freely available, but currently needs to be requested; version 1.1.201 or higher is required   
 - LISFLOOD-FP: version 5.9 extended with BMI-functionality is available at GitLab
 - wflow: the latest version is required for full functionality and can be downloaded from GitHub
 - CaMa-Flood: a BMI'ed version of CaMa-Flood (v3.6.2) is available upon request

## Content of package
 - glofrim: python package containing the functions required to execute the various coupling models
 - scripts: convienence scripts and examples to run GLOFRIM

## Setting up GLOFRIM
We recommend you setup GLOFRIM within its own python environment. You can do so using [conda environments](https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) with the provided envrionment.yml file in the glofrim-py folder. This should also install the required python BMI-wrapper for you.

```
# create environmnet named glofrim
conda env create -f environment.yml
# activate glofrim environment
source activate glofrim
```

To install GLOFRIM do (currently we recommend using the csdms-compliant branch):
```
# get copy of source code from git repos
git clone git@github.com/openearth/glofrim.git@csdms-compliant
# navigate to the py-glofrim folder
cd glofrim/py-glofrim
# install for (the -e links the source code folder for development, this can be left out)
pip install -e <path/to/glofrim>/py-glofrim
```

Then,  download and install the required models. More detailed descriptions how to install the models can be found in the model-specific _bmi.py files. For instance, for PCR-GLOBWB, first install [PC-RASTER](http://pcraster.geo.uu.nl/getting-started/pcraster-on-linux/), then:
```
# pcr-globwb is provided inside the glofrim distribution (will be changed)
pip install <path/to/glofrim>/pcrglobwb_bmi_v203
```

Note that glofrim has only been tested on Linux. 

## Usage
GLOFRIM exists of a series of uniformed BMI wrappers for each model and a overarching BMI wrapper for running coupled models.

To run a coupled model from python use the following lines. The glofrim.ini (see example in root directory) configuration file hold the information of the individual model configuration files and exchanges between the models.
```
from glofrim import Glofrim 
cbmi = Glofrim() # initialize coupled bmi
cbmi.initialize_config(/path/to/glofrim.ini) # initialize the coupling with the glofrim.ini configuration file
```

A basic model run uses the following statements:
```
bmi.get_start_time() # optional: get the model start time
bmi.initialize_model() # initialize model
bmi.update_until(bmi.get_end_time()) # run until set endtime
bmi.finalize()
```

To run stand alone models via the GLOFRIM BMI wrapper you can use followed by the same statements as before:
```
from glofrim import CMF # import the CaMa-Flood bmi wrapper
bmi = CMF(/path/to/model_engine) # intialize bmi with reference to engine (only for non-python models)
bmi.initialize_config(/path/to/model_configuration_file)
```

## Convenience script:
The GLOFRIM library contains a script to run combined and single (for testing purposes) models with a single line from a terminal. This script can be found in the glofirm-py/scripts folder 

GLOFRIM can be executed as follows on Linux command line:
```
python glofrim_runner.py run /path/to/glofrim.ini --env /path/to/glofrim.env -s 200-01-01 -e 2001-01-01
```

For more info on coupled runs, check
```
python glofrim_runner.py run –help
```

and for stand-alone runs:
```
python glofrim_runner.py run_single –help
```

## Literature and sources:
GLOFRIM development and applications \
https://www.geosci-model-dev.net/10/3913/2017/gmd-10-3913-2017.html
		
BMI\
https://bmi-spec.readthedocs.io/en/latest \
http://www.sciencedirect.com/science/article/pii/S0098300412001252

BMI-wrapper for Python\
https://github.com/openearth/bmi-python

PCR-GLOBWB\
https://www.geosci-model-dev.net/11/2429/2018 \
https://github.com/UU-Hydro/PCR_BMI/

Delft3D Flexible Mesh\
https://link.springer.com/article/10.1007%2Fs10236-011-0423-6 \
https://www.deltares.nl/en/software/delft3d-flexible-mesh-suite/#7

LISFLOOD-FP\
http://www.sciencedirect.com/science/article/pii/S002216940000278X \
https://gitlab.com/ChippChapp/LISFLOOD-BMI

CaMa-Flood\
https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2010wr009726 \
https://github.com/hii600/cama-flood_bmi_v3.6.2
	
wflow\
https://wflow.readthedocs.io/en/latest/index.html \
https://github.com/openstreams/wflow

GLOFRIM 1.0 \
https://doi.org/10.5281/zenodo.597107

## Disclaimer:
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

Copyright (C) 2017,2018 Jannis Hoch

The disclaimers/warranty statements of each component involved in this coupling (i.e. PCR-GLOBWB, LIFLOOD-FP, Delft3D Flexible Mesh, BMI Wrapper, CaMa-Flood, wflow)
remain valid unless stated otherwise.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved in case of making use of this coupling script.
And of course don't forget to cite the associated paper as well :)

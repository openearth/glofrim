# GLOFRIM 2.0
Globally Applicable Framework for Integrated Hydrological-Hydrodynamic Modelling (GLOFRIM)

development by Jannis M. Hoch (Utrecht University, Deltares), Dirk Eilander (VU Amsterdam, Deltares), and Hiroaki Ikeuchi (University of Tokyo) \
contact: Jannis M. Hoch (j.m.hoch@uu.nl)

We also want to acknowledge the contributions of all colleagues at the insitues involved in the development of GLOFRIM.

# Manual
GLOFRIM is designed to provide a flexible and modular tool to couple hydrologic, routing, and hydodynamic models across scales.
The coupling process is spatially explicit (i.e. on grid-to-grid basis) and model information is exchanged per (daily) time step.
To establish this coupling scheme, the functions of the Basic Model Interface (BMI) are utilized.

While version 1 allowed for coupling PCR-GLOBWB with either Delft3D Flexible Mesh or LISFLOOD-FP, we extended the framework such that it
now also caters the hydrologic modelling suite WFLOW and the global routing model CaMa-Flood.

With the available models, either 2-step or 3-step step model coupling can be performed.
 - 2-step coupling: hydrology->routing or hydrology->hydrodynamics
 - 3-step coupling: hydrology->routing->hydrodynamics

Currently, the coupling process in only one-directional, i.e. only downstream along the model cascade.
Work is currently performed to extend it to a full feedback loop.

It is important to note that GLOFRIM provides only a tool to coupled models across scales and processes. The quality of simulations therefore still depends on the quality of the model discretizations used.

The framework has successfully been tested on Linux platforms. 
Please note that it running it on Windows is currently not supported.

## Model specs
 - PCR-GLOBWB: since the model does not generically contain BMI function, a bespoke version is provided with the GLOFRIM package.
 - Delft3D Flexible Mesh: the model is freely available, but currently needs to be requested; version 1.1.201 or higher is required   
 - LISFLOOD-FP: version 5.9 extended with BMI-functionality is available at GitLab
 - WFLOW: the latest version is required for full functionality and can be downloaded from GitHub
 - CaMa-Flood: a BMI'ed version of CaMa-Flood (v3.6.2) is available upon request

## Content of package
 - glofrim-py: python package containing the functions required to execute the various coupling models
 - pcrglobwb-bmi_v203: BMI-compatible PCR-GLOBWB model supporting application at 30 arcmin resolution
 - glofrim.ini: example ini-file where the models as well as the exchanged fluxes/paths are specified

## Setting up GLOFRIM
1. Create GLOFRIM python library with ```python setup.py develop``` in glofrim-py
2. Similarly, download and install BMI-wrapper for Python 
3. Download and install models. More detailed descriptions how to install the models can be found in the model-specific _bmi.py files.
4. Run GLOFRIM

## Running the script:
GLOFRIM can be executed as follows on Linux command lines:\
```python glofrim_runner.py run /path/to/glofrim.ini```

Besides, additional arguments can be provided:
 - --env: file containing local paths to model executables
 - -s: start time (yyyy-mm-dd)
 - -e: end time (yyyy-mm-dd)

For more info on coupled runs, check \
```python glofrim_runner.py run –help``` \
and for stand-alone runs: \
```python glofrim_runner.py run_single –help```

## Literature and sources:
GLOFRIM development and applications \
https://www.geosci-model-dev.net/10/3913/2017/gmd-10-3913-2017.html
		
BMI\
https://csdms.colorado.edu/wiki/BMI_Description \
http://www.sciencedirect.com/science/article/pii/S0098300412001252

BMI-wrapper for Python\
https://github.com/openearth/bmi-python

PCR-GLOBWB\
https://www.geosci-model-dev.net/11/2429/2018/

Delft3D Flexible Mesh\
https://link.springer.com/article/10.1007%2Fs10236-011-0423-6 \
https://www.deltares.nl/en/software/delft3d-flexible-mesh-suite/#7

LISFLOOD-FP\
http://www.sciencedirect.com/science/article/pii/S002216940000278X \
https://gitlab.com/ChippChapp/LISFLOOD-BMI

CaMa-Flood\
https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2010wr009726 \
https://github.com/hii600/cama-flood_bmi_v3.6.2
	
WFLOW\
https://wflow.readthedocs.io/en/latest/index.html \
https://github.com/openstreams/wflow

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

The disclaimers/warranty statements of each component involved in this coupling (i.e. PCR-GLOBWB, LIFLOOD-FP, Delft3D Flexible Mesh, BMI Wrapper, CaMa-Flood, WFLOW)
remain valid unless stated otherwise.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved in case of making use of this coupling script.
And of course don't forget to cite the associated paper as well :)

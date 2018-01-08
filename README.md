# GLOFRIM
Globally Applicable Framework for Integrated Hydrological-Hydrodynamic Modelling (GLOFRIM)

developed by Jannis M. Hoch (Utrecht University; j.m.hoch@uu.nl)
(C) 2018

# Content of package
	couplingFramework_v1.py: interface script to execute the coupling process
	default.set: example of the coupling settings-file
	PCRGLOBWB_30min_GLOFRIM_iniFile.ini: example ini-file for PCR-GLOBWB settings to be used with GLOFRIM
	coupling: python library containing functions to perform coupling
	lisflood-bmi-v5.9: BMI-compatible LISFLOOD-FP model at version 5.9
	pcrglobwb-bmi_v203: BMI-compatible PCR-GLOBWB model supporting application at 30 arcmin resolution

# Introduction
The script "couplingFramework_v1.py" can be used to couple PCR-GLOBWB (PCR) with two hydrodynamic models: either Delft3D Flexible Mesh (DFM) or LISFLOOD-FP (LFP).
Please be aware that some variables or settings are filled with placeholders and need to be updated depending on user environment.
The coupling process is facilitated by the use of the Basic Model Interface (BMI).
Therefore it is required that all models to be coupled are (re-)coded such that they are BMI-compatible.
Currently, the coupling process in only one-dimensional, i.e. from PCR to hydrodynamic models.
Work is currently performed to implement a full iterative feedback loop.

So far, the following hydrodynamic models have successfully been coupled and work has been published in scientific journals:
First, Delft3D Flexible Mesh from Deltares which allows for a wide range of hydrodynamic schematizations, e.g. 1D/2D, 2D, flexible mesh or regular grid.
Since DFM is not distributed via GLOFRIM, it needs to be obtained separately. Note, DFM needs to be at version 1.1.201 or higher.
Second, LISFLOOD-FP (version 5.9) from University of Bristol which has been extended with BMI-functionality.
LFP is widely applied in flood hazard and risk assessments and allows for 1D, 2D, and sub-grid channel schematizations, always employing a regular grid.
Please note that the downloadable LISFLOOD-FP version is not meant for further unauthorized distribution.

The "default.set" file provided in the folder is a template where all required paths to the models (so-files for Linux and dll-files for windows) have to be provided.
Besides, all model settings can be set there.
The framework has successfully been tested on Linux platforms. Please note that running GLOFRIM on Windows is currently only supported for DFM as no dll is compiled for LFP.

Since the ini-file of PCR to be used with GLOFRIM differs from the default, a GLOFRIM-compatible template is provided as well.
All PCR related settings can be defined here. Further information can be found in the file itself.

Please note:
While it's the purpose that the couplingFramework_v1.py file is generic, a separate set-file, ini-file, and mdu/par-file are required for each GLOFRIM set-up and run.
In the set-file, particularly the paths to the ini-file and the mdu/par-file need to be set depending on model schematizations chosen.
In the ini-file, the output folder for PCR output needs to be different per run, at least if multiple runs are ought to be executed simultaneously.
In the mdu/par-file, the output folder for DFM/LFP output needs to be different per run, at least if multiple runs are ought to be executed simultaneously.

# Getting started
Before you can apply GLOFRIM, a few preparatory steps are necessary.
First, the packages "coupling" and "pcrglobwb-bmi_v203" need to be converted to python libraries by executing "python setup.py develop" in the respective folders.
You can check whether they are compiled correctly by importing "coupling_PCR_DFM" and "pcrglobwb_bmi_v203", respectively, into your python environment.
Besides, a python-compatible BMI-wrapper needs to be downloaded (see link below) and also converted to a python library.
Again, you can test the correct installation by importing "bmi.wrapper" into your python environment.
In case PCR shall be coupled to DFM, a DFM version has to be installed first. Currently please contact Deltares (see link below) to receive one. Please specify which environment you are
working it (Linux or Windows) and also state your purpose.
After having received the model, please set the path to the dflowfm.dll-file (Windows) or the libdflowfm.so-file (Linux) in couplingFramework_v1.py.
In case PCR shall be coupled to LFP, please cd to the folder lisflood-bmi-v5.9 and compile the programm my executing "make".
In case you encounter any issues, please consider adjusting the makefile according to your Linux environment.
Then set the path to liblisflood.so (Linux) in couplingFramework_v1.py.
Ideally, the paths to the models remain unaltred which allows for using the couplingFramework_v1.py for all GLOFRIM runs.

For questions, lessons learnt, experiences made or if any problems are encountered, please don't hesitate to contact me (j.m.hoch@uu.nl)

# Literature and sources:
BMI
Conceptual description: https://csdms.colorado.edu/wiki/BMI_Description
Peckham et al, 2013: http://www.sciencedirect.com/science/article/pii/S0098300412001252
Download page BMI-wrapper: https://github.com/openearth/bmi-python

PCR-GLOBWB
Technical details: http://vanbeek.geo.uu.nl/suppinfo/vanbeekbierkens2009.pdf
Sutanudjaja et al, 2017: https://www.geosci-model-dev-discuss.net/gmd-2017-288/

Delft3D Flexible Mesh
Technical description: https://link.springer.com/article/10.1007%2Fs10236-011-0423-6
Contact page: https://www.deltares.nl/en/software/delft3d-flexible-mesh-suite/#7

LISFLOOD-FP
Bates et al, 2010: http://www.sciencedirect.com/science/article/pii/S002216940000278X
Bates et al, 2010: http://www.sciencedirect.com/science/article/pii/S0022169410001538

APPLICATIONS OF GLOFRIM
Hoch et al, 2017: https://www.hydrol-earth-syst-sci.net/21/117/2017/hess-21-117-2017.html
Hoch et al, 2017: https://www.geosci-model-dev.net/10/3913/2017/

# Running GLOFRIM:
To run the script, an set-file containing the required specifications and paths is necessary.
Using python, run this file along with the set-file as follows:
	python couplingFramework_v1.py default.set

# Disclaimer:
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

Copyright (C) 2018 Jannis Hoch

The disclaimers of each component involved in this coupling (i.e. PCR-GLOBWB, LIFLOOD-FP, Delft3D Flexible Mesh, BMI Wrapper)
remain valid unless stated otherwise.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved in case of making use of this coupling script.
And of course don't forget to cite the associated paper as well :)

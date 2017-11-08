# GLOFRIM
Globally Applicable Framework for Integrated Hydrological-Hydrodynamic Modelling (GLOFRIM)

# NESTING
This branch will be used for developed of a coupled hydrology -> routing -> 1d/2d hydrodynamics model framework. 

# Content of package
	couplingFramework_v1.py: script to execute the coupling process
	default.set: example of the coupling settings-file
	PCRGLOBWB_30min_GLOFRIM_iniFile.ini: example file to be used with GLOFRIM for PCR-GLOBWB settings
	coupling: python library containing functions to perform coupling
	lisflood-bmi-v5.9: BMI-compatible LISFLOOD-FP model at version 5.9
	pcrglobwb-bmi_v203: BMI-compatible PCR-GLOBWB model supporting application at 30 arcmin resolution

# Manual
The script provided can be used to couple PCR-GLOBWB with two hydrodynamic models: Delft3D Flexible Mesh (DFM) and LISFLOOD-FP (LFP).
The coupling process is facilitated by the use of the Basic Model Interface (BMI).
Currently, the coupling process in only one-dimensional, i.e. from PCR-GLOBWB to hydrodynamic models.
Work is currently performed to extend it to a full feedback loop.

So far, the following hydrodynamic models have successfully been coupled:
First, the in-house Delft3D Flexible Mesh which needs to be at version 1.1.201 or higher.
Second, the LISFLOOD-FP model from University of Bristol at version 5.9 which has been extended with BMI-functionality.
Please note that the downloadable LISFLOOD-FP version is not meant for further unauthorized distribution.

The set-file provided in the folder is a template where all required paths needs to be set.
All relevant information regarding model set-ups and settings are provided there.
The framework has successfully been tested on Linux platforms. Please note that it running it on Windows is currently not supported as no dll is yet compiled for LFP.

# Getting started
Before coupling is possible, a few steps need to be taken.
First, the packages "coupling" and "pcrglobwb-bmi_v203" need to be converted
to python library by typing "python setup.py develop" in the respective folders. Besides, a python-compatible BMI-wrapper needs to be downloaded (see link below) and also converted
to a python library.
In case PCR shall be coupled to DFM, a DFM version has to be installed first. Currently please contact Deltares (see link below) for receiving one. Please specify which environment you are
working it (Linux or Windows) and also state your purpose.
After having received the data, please set the path to the dflowfm.dll-file (Windows) or the libdflowfm.so-file (Linux) in couplingFramework_v1.py.
In case PCR shall be coupled to LFP, please cd to the folder lisflood-bmi-v5 and "make" the programme. In case you encounter any issues, please consider adjusting the makefile.
Then set the path to liblisflood.so (Linux) in couplingFramework_v1.py.

For questions, lessons learnt, experiences made or if any problems are encountered, please contact Jannis Hoch (j.m.hoch@uu.nl)

# Literature and sources:
	BMI
	https://csdms.colorado.edu/wiki/BMI_Description
	http://www.sciencedirect.com/science/article/pii/S0098300412001252

	BMI-wrapper
	https://github.com/openearth/bmi-python

	PCR-GLOBWB
	http://vanbeek.geo.uu.nl/suppinfo/vanbeekbierkens2009.pdf

	Delft3D Flexible Mseh
	https://link.springer.com/article/10.1007%2Fs10236-011-0423-6 (Technical Description)
	https://www.deltares.nl/en/software/delft3d-flexible-mesh-suite/#7 (Contact page)

	LISFLOOD-FP
	http://www.sciencedirect.com/science/article/pii/S002216940000278X

# Running the script:
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

Copyright (C) 2017 Jannis Hoch

The disclaimers of each component involved in this coupling (i.e. PCR-GLOBWB, LIFLOOD-FP, Delft3D Flexible Mesh, BMI Wrapper)
remain valid unless stated otherwise.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved in case of making use of this coupling script.
And of course don't forget to cite the associated paper as well :)

# GLOFRIM
Globally Applicable Framework for Integrated Hydrological-Hydrodynamic Modelling (GLOFRIM)

# Content of package
The downloadable release contains the following files:
	couplingFramework_v1.py: script to execute the coupling process
	default.ini: example of the coupling settings-file
	coupling: python library containing functions to perform coupling
	lisflood-bmi-v5.9: BMI-compatible LISFLOOD-FP model at version 5.9
	pcr-globwb-203-30min-1way-prefactored: BMI-compatible PCR-GLOBWB model supporting application at 30 arcmin resolution
	
# Manual
The script provided can be used to couple PCR-GLOBWB with two hydrodynamic models: Delft3D Flexible Mesh (DFM) and LISFLOOD-FP (LFP).
The coupling process is facilitated by the use of the Basic Model Interface (BMI). 
Currently, the coupling process in only one-dimensional, i.e. from PCR-GLOBWB to hydrodynamic models. 
Work is currently performed to extend it to a full feedback loop.

So far, the following hydrodynamic models have successfully been coupled:
First, the in-house Delft3D Flexible Mesh which needs to be at version 1.1.201 or higher.
Second, the LISFLOOD-FP model from University of Bristol at version 5.9 which has been extended with BMI-functionality.
Please note that the downloadable LISFLOOD-FP version is not meant for further unauthorized distribution.

The ini-file provided in the folder is a template where all required paths needs to be set.
All relevant information regarding model set-ups and settings are provided there.
The framework has successfully been tested on Linux platforms. Please note that it running it on Windows is limited as no dll is yet compiled for LFP.

# Getting started
Before coupling is possible, a few steps need to be taken. 
First, the packages "coupling" and "pcr-globwb-203-30min-1way-prefactored" need to be converted
to python library by typing "python setup.py develop" in the respective folders. Besides, the BMI-wrapper needs to be downloaded (see link below) and also converted
to a python library.
Additionally, a DFM version has to be installed. Currently please contact Deltares (see link below) for receiving one. Please specify which environment you are
working it (Linux or Windows) and also state your purpose.

For questions, lessons learnt, experiences made or if any problems are encountered, please contact Jannis Hoch (j.m.hoch@uu.nl)

# Literature and sources:
	BMI
	https://csdms.colorado.edu/wiki/BMI_Description
	http://www.sciencedirect.com/science/article/pii/S0098300412001252
	
	bmi.wrapper
	https://github.com/openearth/bmi-python
	
	PCR-GLOBWB
	http://vanbeek.geo.uu.nl/suppinfo/vanbeekbierkens2009.pdf
	
	Delft3D Flexible Mseh
	https://link.springer.com/article/10.1007%2Fs10236-011-0423-6 (Technical Description)
	https://www.deltares.nl/en/software/delft3d-flexible-mesh-suite/#7 (Contact page)
	
	LISFLOOD-FP
	http://www.sciencedirect.com/science/article/pii/S002216940000278X

# Running the script:
To run the script, an ini-file containing the required specifications and paths is necessary.
Using python, run this file along with the ini-file as follows:
	python couplingFramework_v1.py default.ini
	
# Disclaimer:
The disclaimers of each component involved in this coupling (i.e. PCR-GLOBWB, LIFLOOD-FP, Delft3D Flexible Mesh, BMI Wrapper)
remain valid unless stated otherwise.
No warranty/responsibility for any outcome of using this coupling script.
Please ensure to cite the models involved in case of making use of this coupling script.
And of course don't forget to cite the associated paper as well :)


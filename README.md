# GLOFRIM
Globally Applicable Framework for Integrated Hydrological-Hydrodynamic Modelling (GLOFRIM)

# MANUAL COUPLING FRAMEWORK v1.0
The files in this folder can be used to couple PCR-GLOBWB with two hydrodynamic models: Delft3D Flexible Mesh and LISFLOOD-FP.
The coupling process is facilitated by the use of the Basic Model Interface (BMI). 
Currently, the coupling process in only one-dimensional, i.e. from PCR-GLOBWB to hydrodynamic models. 
Work is currently performed to extend it to a full feedback loop.

So far, the following hydrodynamic models have successfully been coupled:
First, the in-house Delft3D Flexible Mesh (FM) which needs to be at version 1.1.201 or higher.
Second, the LISFLOOD-FP model from University of Bristol at version 5.9 which has been extended with BMI-functionality.
Please note that the downloadable LISFLOOD-FP version is not meant for further unauthorized distribution! 

The ini-file provided in the folder is a template where all required paths needs to be set.
All relevant information regarding model set-ups and settings are provided there.
The framework has successfully been tested on Linux platforms. Please note that it running it on Windows is limited as no dll is yet compiled for LFP.

Together with the code for the computational framework and the settings-file, we provide the code for LISFLOOD-FP 5.9 including BMI and PCR-GLOBWB with BMI.
Please be aware that due to on-going model development these versions do not represent the latest versions.
Yet, it is ensured that all key model functionalities are accessible.

Along with the files required for coupling, we provide explanatory model files for the three models currently implemented.
All files should be self-explanatory. In case issues with model set-up are encountered, we refer to the specific read-mes and manuals of the models. 

For questions, lessons learnt, experiences made or if any problems are encountered, please contact Jannis Hoch (j.m.hoch@uu.nl)

# Literature and sources:
	BMI         -> https://csdms.colorado.edu/wiki/BMI_Description
				-> http://www.sciencedirect.com/science/article/pii/S0098300412001252
	bmi.wrapper -> https://github.com/openearth/bmi-python
	PCR-GLOBWB	-> http://vanbeek.geo.uu.nl/suppinfo/vanbeekbierkens2009.pdf
	Delft3D FM	-> https://link.springer.com/article/10.1007%2Fs10236-011-0423-6
	LISFLOOD-FP	-> http://www.sciencedirect.com/science/article/pii/S002216940000278X

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


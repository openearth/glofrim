# GLOFRIM
Globally Applicable Framework for Integrated Hydrological-Hydrodynamic Modelling (GLOFRIM)

# Introduction:
The coupling is achieved by making use of the Basic Model Interface (BMI) which allows for initializing,
updating, data manipulation, and finalization of models from a shell-like environment as for instance this python-script. 
For couple Delft3D FM or LISFLOOD-FP, the python-module "bmi.wrapper" has to be loaded.
For Delft3D FM, any compiled version (>1.1.201.48898) has already implemented a BMI-compatible structure, and the 
required variables accessible to the user (downloadable upon request).
For LISFLOOD-FP, however, a specifically designed version needs to be compiled which is currently only available for
version 5.9 as the model is not generically BMI-compatible (also downloadable here).
Also for PCR-GLOBWB, a BMI-compatible version needs to be used (also downloadable here).

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


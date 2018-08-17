.. requirements:

*******************************
GLOFRIM requirements and set-up
*******************************

The tool was entirey written in Python 2.7 and tested on Linux systems.

.. note::
    Migration to Python 3.x is envisaged but not scheduled yet.

.. note::
    GLOFRIM requires `PCRaster <http://pcraster.geo.uu.nl/getting-started/pcraster-on-linux/>`_ to be installed on your machine!
    Please ensure this is the case before proceeding with the GLOFRIM set-up.

.. warning::
    GLOFRIM is not supported on Windows systems!

To install and run GLOFRIM on your Linux machine, several python packages are required.
For convenience, we recommend to set up GLOFRIM within its own python environment.
You can do so using `conda environments <https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ 
with the provided envrionment.yml file in the glofrim-py folder.

This should also install the required python BMI-wrapper for you::

    # create environmnet named glofrim
    conda env create -f environment.yml
    # activate glofrim environment
    source activate glofrim    

To install GLOFRIM do::

    # get copy of source code from git repos
    git clone git@github.com/openearth/glofrim.git@csdms-compliant
    # navigate to the py-glofrim folder
    cd glofrim/py-glofrim
    # install for (the -e links the source code folder for development, this can be left out)
    pip install -e <path/to/glofrim>/py-glofrim

.. note::
    The currently most stable branch is csdms-compliant. As development and merging goes along, this may change however.
    Any change will be documented here.

Once the glofrim environment is successfully installed, the actual models can be installed.




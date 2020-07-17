.. _requirements:

*******************************
GLOFRIM requirements and set-up
*******************************

.. note::
    GLOFRIM is available for Python 3.x on GitHub. The v2.0 release on Zenodo is still based on Python 2.7 though!

To install and run GLOFRIM on your Linux machine, several python packages are required.
For convenience, we recommend to set up GLOFRIM within its own python environment.
You can do so using `conda environments <https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ 
with the provided envrionment.yml.

This should also install the required python BMI-wrapper for you.

.. code-block:: console

    conda env create -f environment.yml
    conda activate glofrim    

To install GLOFRIM do:

.. code-block:: console

    git clone git@github.com/openearth/glofrim.git
    cd glofrim/py-glofrim
    pip install -e <path/to/glofrim>/py-glofrim

Once the glofrim environment is successfully set up, the actual models can be installed.

.. note::

    GLOFRIM installation via pip is currently not yet supported.

Note for installing PCR and WFL
-------------------------------
Both hydrological models use the `PCRaster <http://pcraster.geo.uu.nl/getting-started/pcraster-on-linux/installation-linux/>`_ library within its modelling framework. 
For a smooth GLOFRIM experience, please ensure PCRaster has been installed before installing PCR or WFL. By default this happens while setting up the environment, but better double-check.

WFL is also directly installed within the environment. PCR, however, must be downloaded and installed separately.
    
Note for installing CMF, LFP and DFM
------------------------------------
All hydrodynamic models are written in other languages than Python, such as Fortran anc C++. These models need to be compiled on your linux machines, check the
documentation of the individual models for more info.

The paths to the share libraries of each need to be provided to GLOFRIM. These can be set in the GLOFRIM ini file or provided in a file called environment.env 
in the glofrim root dir in the *engines* section, see :ref:`The GLOFRIM configuration file`.

.. note::

    In addition to the source code, to set up local models and to generate the interpolation matrix for runoff inputs to CMF, additional Fortran scripts need to be compiled.
    GLOFRIM requires the generate_inpmat.F script to setup the interpolation scheme based on domain of the upstream model in the model cascade. 
    Carefully follow the steps laid out in the README.txt in the glofim/script/generate_inpmat folder for more info.

Testing your installation
-------------------------
Once you've followed all steps above: installed GLOFRIM; the models you want to work with and you have set the environment.env file,
you can test your setup by running the unit tests provided in glofrim/tests folder. 

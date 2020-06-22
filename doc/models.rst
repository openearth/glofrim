.. models:

******************
Model descriptions
******************

Supported models in GLOFRIM
===========================

In the current version, five different models are supported by GLOFRIM: two hydrologic models (PCR-GLOBWB and WFLOW), one routing model (CaMa-Flood), and
two hydrodynamic models (Delft3D Flexible Mesh and LISFLOOD-FP).

PCR-GLOBWB
----------
PCR-GLOBWB (hereafter PCR) is a global hydrology and water resources model, fully integrating water use. Sector-specific water demand, groundwater and surface water withdrawal, 
water consumption, and return flows are dynamically calculated at every time step and interact directly with the simulated hydrology.

The model simulates vertical exchanges between surface, two soil moistere layers, and a groundwater layer. Exchange precipitation is, after accounting for irrigatoin and water use abstractions,
routed along a local drainage network using the kinematic wave approximation.

For more information about PCR, we refer to [Sutanudjaja2018]_.

Model availability
^^^^^^^^^^^^^^^^^^^

The default and published stable version of PCR is available at `Zenodo <https://doi.org/10.5281/zenodo.595656>`_ which is, however, not suitable for application within GLOFRIM.

To apply PCR within GLOFRIM, a separate model version of PCR-GLOBWB needs to be used. It can be found `here <https://doi.org/10.5281/zenodo.3345900>`_.

.. note::

    Due to the evolving nature of model code and thus also of PCR-GLOBWB, recently added functionality of 
    PCR-GLOBWB may not be available in the here downloadable version.
    If you encounter any issues, please contact the developers!

.. note::

    PCR within GLOFRIM was successfully applied at 30 arc-min spatial resolutions. 
    The 05 arc-min version or schematization at any finer spatial resolution were so far only used in test model and not for detailed studies.

WFLOW
-----
WFLOW (hereafter WFL) is a framework for distributed hydrologic modelling, currently including the following hydrologic models: SBM, HBV, GR4, W3RA and PCRGLOB-WB.
In addition, routing can be simulated using the kinematic wave approximation.
Depending on the chosen model, the processes described differ and therefore we do not provide a more detailed description of the model routines.

Despite the different models available, WFL ensures that the main model properties and structure (nomenclature, grid properties etc.) are mantained.

For a full documentation of WFL code, models, application, and other functionalities, please visit the `WFLOW readthedocs <https://wflow.readthedocs.io/en/latest/>`_ page.

Model availability
^^^^^^^^^^^^^^^^^^^

WFL contains a BMI adapter by default and therefore no alternative version has to be applied to use it within GLOFRIM.

All code is available at `WFLOW GitHub <https://github.com/openstreams/wflow/>`_ and distributed under the GNU GPL 3.0 license.

.. note::

    For GLOFRIM, only the SBM and W3RA models were tested successfully. Application of other WFL models should be straightforward, but that is uncharted territory at the moment.

CaMa-Flood
----------
CaMa-Flood (hereafter CMF) is designed to simulate the hydrodynamics in continental-scale rivers. It has global extent and routes water along river networks solving the 
local inertia equation. Water level and flooded area are determined from the water storage at each unit catchment. Since water balance per unit catchment is the only prognostic 
variable, simulation are efficient.

Since CMF is a global model, only little additional work has to be done to obtain smaller discretizations. Also, this faciliates comparability between test studies since differences
in input data are marginal.

For further model description and application, please see [Yamazaki2011]_ as well as the model's `manual <http://hydro.iis.u-tokyo.ac.jp/~yamadai/cama-flood/Manual_CaMa-Flood_v362.pdf>`_.

.. note::

    For the development of GLOFRIM, we applied version 3.6.2. Higher version may differ.

Model availability
^^^^^^^^^^^^^^^^^^^

The original non bmi'ed model code is available upon request `here <http://hydro.iis.u-tokyo.ac.jp/~yamadai/cama-flood/>`_.

.. note::

    To obtain a bmi'ed version (v.3.6.2) to be applied within GLOFRIM, please contact the developers.

Delft3D Flexible Mesh
---------------------
Delft3D Flexible Mesh (hereafter DFM) is a fully fledged hydrodynamic model, solving the full shallow water equations. Models can be set up in 1D, 2D, 3D or 1D/2D. For large-scale
inundation modelling, a 1D/2D set-up is preferrable and thus also the default design supported by GLOFRIM.

While DFM is used for inundation modelling within GLOFRIM, it also supports other applications, such as salt intrusion or sediment transport simulations. In principle, the
coupled simulations by GLOFRIM could be extended with other morphogeohydrologic variables, further increasing the representation of relevant processes.

As the name indicates, DFM allows for discretizing the model domain with a flexible mesh; that is, varying cell geometry and size. By using flexible meshes, the representation
of topographical features such as as river bends can be realized with fine spatial resolution, while areas with less dynamic processes require only coarser grids. As a result,
the run time needed can be reduced significantly.

Since DFM genercially contains a BMI adapter, any DFM version can be used as long as it satisfies the minimal version requirement.

For documentation of the model scheme, [Kernkamp2011]_ provides the theoretical and computatonal background. Besides, extensive `technical <https://content.oss.deltares.nl/delft3d/manuals/D-Flow_FM_Technical_Reference_Manual.pdf>`_ 
and `user <https://content.oss.deltares.nl/delft3d/manuals/D-Flow_FM_User_Manual.pdf>`_ manuals are available.

Model availability
^^^^^^^^^^^^^^^^^^^

The model is free to use, but currently not yet openly available. The DFM development team needs to be contacted for a DFM version.
Please see the `DFM website <https://oss.deltares.nl/web/delft3dfm/home>`_ for contact information.

.. note::

    DFM version higher than 1.1.201 is required to work with GLOFRIM, the framework has successfully been tested with version 1.1.201.

LISFLOOD-FP
-----------
LISFLOOD-FP (hereafter LFP) is a well tested and widely used hydrodynamic model specifically designed to simulate floodplain inundation in a computationally efficient manner over complex topography. 
It computes water depths in each grid cell at each time step, and hence can simulate the dynamic propagation of flood waves over fluvial, coastal, and estuarine floodplains.

While LFP also allows for 1D and 2D set-ups, only the sub-grid channel design was employed within GLOFRIM due to is improved accuracy.

A major advantage of LFP is its easy model creation which requires, for the simplest set-up, only ascii files describing the DEM, the channel width and bed level elevation, as
well as the river bank height. The computational grid is regular in all applications.

The initial paper documenting LFP's computational scheme is [Bates2010]_. More model and background information can be found on the `LISFLOOD-FP <https://www.bristol.ac.uk/geography/research/hydrology/models/lisflood/>`_ website.

Model availability
^^^^^^^^^^^^^^^^^^^

The bmi'ed version of LFP (v. 5.9) can freely be downloaded from `Zenodo <https://doi.org/10.5281/zenodo.1479836>`_. 
A test version of the default model can be requested via this `form <https://www.bristol.ac.uk/geography/research/hydrology/models/lisflood/downloads/>`_.

.. note::

    The downloadable bmi'ed version is based on LFP version 5.9 and not updated with recent updates.
    The computational scheme is, nevertheless, identical and inundation simulations are not affected.

Adding new models
-----------------
It's (relatively) easy to extend GLOFRIM with new models.
A requirement is that the model to be added contains BMI functions and follows the conventions used in the python-BMI files of the other models.
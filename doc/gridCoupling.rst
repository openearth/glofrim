.. _gridCoupling:

*******************************
Coupling different models
*******************************
When different models with different scopes are coupled, this may in a range of possible coupling
designs. To cater for this flexibility, GLOFRIM allows for various ways of coupling grids and exchanging
model information.

To establish online and spatially explicit coupling between models, we employ the `Basic Model Interface <https://csdms.colorado.edu/wiki/BMI_Description>`_.
An elaborated outline of this concept is provided at :ref:`basicModelInterface`.

Spatial coupling of model grids
===============================

.. _2dto2d:

2D to 2D
--------
For models containing only of 2D model grid, meshes, or unit catchments, for instance PCR and CMF.

By means of the developed grid API, GLOFRIM detects the corresponding cells of different model grids and
assigns a 1-n index list which is subsequently used to attribute the right volumes to the right grids.

If n>1, then the volumes from the providing model are devided equally over all coupled cells in the receiving model
to maintain a correct water balance between models.

2D to 1D
--------
If coupling a 2D grid/mesh/unit catchment to the channel network of either DFM or LFP.

The workflow is identical as for :ref:`2dto2d`, but before the 1-n index list is created, GLOFRIM finds those indices in the 
receiving model that belong to river channels. To that end, GLOFRIM employs model-specific properties which are unique to river
channels and can therefore be used for the separation.

.. note::
    GLOFRIM aims at using the default input files of the models to perform the spatial coupling of grids and networks. For
    CMF, however, it is still needed to convert bin/ctl-files to geoTiff-files, particularly the nextxy.bin file for look-up of
    the flow direction.
    We provide a script to convert bin/ctl-files to tif-files.::
        python <path/to>/ cama_maps_io.py <map/folder>/nextxy.ctl. 

Exchange of model information
=============================

.. todo::
    Insert figures for illustration

With GLOFRIM, model variables can be exchanged for the entire grid/network of a model or only for the most upstream
cells/nodes, depending on model set-up and envisaged application.

To define the way model information is exchanged, this has to be specified in the ini-file as explained in :ref:`the_ini_file`.

.. note::
    The way model information is exchanged may differ per variable!

Entire domain
-------------
If the extent of the models domains of the different models match, information can be exchanged over the entire domain.
Alternatively, this option may be chosen in a nested setting if states are exchanged that were not routed in the upstream model.

Only upstream
-------------
In case the downstream model is nested into the upstream model (i.e. the extent is smaller), the upstream model
may already route discharge until the edge of the downstream model is reached.
This may happen if a hydrodynamic model only represents the river delta while the routing and/or hydrologic model capture
the entire basin.
To account for the processes already simulated upstream, the resulting fluxes/states at the edge must only there be added
to the downstream model.






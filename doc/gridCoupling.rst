.. _gridCoupling:

*******************************
Coupling different models
*******************************
When different models with different scopes are coupled, this may in a range of possible coupling
designs. To cater for this flexibility, GLOFRIM allows for various ways of coupling grids and exchanging
model information.

Spatial coupling of model grids
===============================

2D to 2D
--------
For models containing only of 2D model grid, meshes, or unit catchments, for instance PCR and CMF.

2D to 1D
--------
If coupling a 2D grid/mesh/unit catchment to the channel network of either DFM or LFP.

.. todo::
    add more technical detail on how grids are coupled.

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






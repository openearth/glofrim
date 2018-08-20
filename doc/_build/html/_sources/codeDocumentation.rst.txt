.. codeDocumentation:


**************************
GLOFRIM code documentation
**************************

Initialization
==============
Before a coupled model can be run, the individual model have to be initialized and all exchanges
between them must be defined.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.initialize_config
   Glofrim.initialize_model
   Glofrim.initialize
   Glofrim.set_exchanges
   Glofrim.set_out_dir

Execution
===============
GLOFRIM provides a range of functions to retrieve and change time information,
update models, as well as exchange content between models during execution.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.get_start_time
   Glofrim.get_current_time
   Glofrim.get_time_step
   Glofrim.update
   Glofrim.update_until
   Glofrim.exchange
   Glofrim.exchange_same_grid
   Glofrim.exchange_at_indices

Variable Getter and Setter Functions
------------------------------------
Functions to retrieve and overwrite values for either entire grid or 
only certain indices.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.get_value
   Glofrim.get_value_at_indices
   Glofrim.set_value
   Glofrim.set_value_at_indices

Finalization
============
It is possible to change model end times and finalize model states after execution.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.get_end_time
   Glofrim.finalize

Auxiliary functions
===================
In addition, there are several auxiliary functions built in GLOFRIM to check states and properties
of models as well as their components, variables, and attributs.

Model Information Functions
---------------------------
General information about model structure and variables.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.get_model_type
   Glofrim.get_component_name
   Glofrim.get_input_var_names
   Glofrim.get_output_var_names

Variable Information Functions
------------------------------
Providing information about properties 
of model variables.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.get_var_type
   Glofrim.get_var_units
   Glofrim.get_var_rank
   Glofrim.get_var_size
   Glofrim.get_var_shape
   Glofrim.get_var_nbytes
   Glofrim.get_time_units

Grid Information Functions
---------------------------
Information about model grid.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.get_grid_type

Attribute/ Config Information Functions
---------------------------------------
Functions providing insights in settings of model 
configuration files.

.. automodule:: glofrim

.. currentmodule:: glofrim

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Glofrim.set_start_time
   Glofrim.set_end_time
   Glofrim.get_attribute_names
   Glofrim.get_attribute_value
   Glofrim.set_attribute_value
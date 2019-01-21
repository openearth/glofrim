#! /usr/bin/env python

"""
To create this file the original version of the CSDMS BMI Python Language Binding
from https://github.com/csdms/bmi-python/blob/master/bmi/bmi.py was extended with a number of extra functions
so that it can be used in OpenDA. All extra functions are contained in the class EBmi.
Also added additional comments. Where the original version of the CSDMS BMI Python Language Binding was ambiguous, the
information from http://csdms.colorado.edu/wiki/BMI_Description and common sense were used to fill in most of the gaps.


https://github.com/eWaterCycle/bmi/blob/master/src/main/python/bmi.py
"""

from abc import ABCMeta, abstractmethod

class Bmi(object):
    """
    Interface (abstract base class) for a model that implements the CSDMS BMI (Basic Model Interface).
    """

    __metaclass__ = ABCMeta
    grid = None
    _name = 'dummy'
    _long_name = 'dummy'
    _var_units = {'var1': 'm/s', 'var2': 'm'}
    _input_var_names = ['var1']
    _output_var_names = ['var2']

    """
    Model Control Functions
    """

    @abstractmethod
    def initialize(self, filename):
        """
        Initialize the model.
        Input parameters:
        File filename: path and name of the configuration file for the model.
        """
        pass

    @abstractmethod
    def update(self):
        """
        Update the model to the next time step.
        """
        pass

    @abstractmethod
    def update_until(self, time):
        """
        Update the model until the given time.
        Input parameters:
        double time: time in the units and epoch returned by the function get_time_units.
        """
        pass

    # @abstractmethod
    # def update_frac(self, time_frac):
    #     """
    #     ???
    #     Input parameters:
    #     double time_frac: ???
    #     """
    #     pass

    @abstractmethod
    def finalize(self):
        """
        Finalize the model.
        """
        pass

    """
    Model Information Functions
    """

    def get_component_name(self):
        """
        Return value:
        String: identifier of the model.
        """
        return self._name

    def get_input_var_names(self):
        """
        Return value:
        List of String objects: identifiers of all input variables of the model.
        """
        return self._input_var_names

    def get_output_var_names(self):
        """
        Return value:
        List of String objects: identifiers of all output variables of the model.
        """
        return self._output_var_names

    """
    Variable Information Functions
    """


    def get_var_type(self, long_var_name):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Return value:
        String: data type of the values of the given variable, e.g. Numpy datatype string.
        """
        return str(self.get_value(long_var_name).dtype)

    def get_var_units(self, long_var_name):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Return value:
        String: unit of the values of the given variable. Return a string formatted using the UDUNITS standard from Unidata.
        """
        return self._var_units.get(long_var_name, 'unknown')

    def get_var_rank(self, long_var_name):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Return value:
        Integer: number of dimensions of the given variable.
        """
        return self.get_value(long_var_name).ndim

    def get_var_size(self, long_var_name):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Return value:
        Integer: total number of values contained in the given variable, e.g. gridCellCount.
        """
        return self.get_value(long_var_name).size

    def get_var_nbytes(self, long_var_name):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Return value:
        ???: ???
        """
        return self.get_value(long_var_name).nbytes

    @abstractmethod
    def get_start_time(self):
        """
        Return value:
        double: start time of the model in the units and epoch returned by the function get_time_units.
        """
        pass

    @abstractmethod
    def get_current_time(self):
        """
        Return value:
        double: current time of the model in the units and epoch returned by the function get_time_units.
        """
        pass

    @abstractmethod
    def get_end_time(self):
        """
        Return value:
        double: end time of the model in the units and epoch returned by the function get_time_units.
        """
        pass

    @abstractmethod
    def get_time_step(self):
        """
        Return value:
        #NOTE: changed from double type to datetime.timedelta type for convenience
        timedelta: duration of one time step as timedelta object
        """
        pass

    @abstractmethod
    def get_time_units(self):
        """
        Return value:
        String: unit and epoch of time in the model. Return a string formatted using the UDUNITS standard from Unidata.
        """
        pass

    """
    Variable Getter and Setter Functions
    """

    @abstractmethod
    def get_value(self, long_var_name):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Return value:
        Numpy array of values in the data type returned by the function get_var_type: all values of the given variable.
                                                                                      For a 2D grid these values must be in row major order, starting with the bottom row.
        """
        pass

    @abstractmethod
    def get_value_at_indices(self, long_var_name, inds):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        List of Lists of integers inds: each nested List contains one index for each dimension of the given variable,
                                        i.e. each nested List indicates one element in the multi-dimensional variable array,
                                        e.g. [[0, 0, 0], [0, 0, 1], [0, 15, 19], [0, 15, 20], [0, 15, 21]] indicates 5 elements in a 3D grid.
                                        For a grid the indices start counting at the grid origin, i.e. the lower left corner for a 2D grid.
        Return value:
        Numpy array of values in the data type returned by the function get_var_type: one value for each of the indicated elements.
        """
        pass

    @abstractmethod
    def set_value(self, long_var_name, src):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Numpy array of values in the data type returned by the function get_var_type src: all values to set for the given variable.
                                                                                          For a 2D grid these values must be in row major order, starting with the bottom row.
        """
        pass

    @abstractmethod
    def set_value_at_indices(self, long_var_name, inds, src):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        List of Lists of integers inds: each nested List contains one index for each dimension of the given variable,
                                        i.e. each nested List indicates one element in the multi-dimensional variable array,
                                        e.g. [[0, 0], [0, 1], [15, 19], [15, 20], [15, 21]] indicates 5 elements in a 2D grid.
                                        For a grid the indices start counting at the grid origin, i.e. the lower left corner for a 2D grid.
        Numpy array of values in the data type returned by the function get_var_type src: one value to set for each of the indicated elements.
        """
        pass

    
    # NOTE all methods are replaced by Grid object methods using rasterio, pyugrid etc., see grid.py
    # """
    # Grid Information Functions
    # """
    # @abstractmethod
    # def get_grid.type(self, long_var_name):
    #     """
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     BmiGridType type of the grid geometry of the given variable.
    #     """
    #     pass
    # @abstractmethod
    # def get_grid_shape(self, long_var_name):
    #     """
    #     Return value:
    #     tuple of integers: the sizes of the dimensions of the given variable, e.g. [400, 500] for a 2D grid with 400 rows and 500 columns.
    #                       The dimensions are ordered [y, x] or [z, y, x].
    #     """
    #     pass

    # @abstractmethod
    # def get_grid_spacing(self, long_var_name):
    #     """
    #     Only return something for variables with a uniform grid. Otherwise raise ValueError.
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     List of doubles: the size of a grid cell for each of the dimensions of the given variable, e.g. [cellHeight, cellWidth] for a 2D grid.
    #                      The dimensions are ordered [y, x] or [z, y, x].
    #     """
    #     pass

    # @abstractmethod
    # def get_grid_origin(self, long_var_name):
    #     """
    #     Only return something for variables with a uniform grid. Otherwise raise ValueError.
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     List of doubles: the coordinate of the grid origin for each of the dimensions of the given variable.
    #                      For a 2D grid this must be the lower left corner of the grid.
    #                      The dimensions are ordered [y, x] or [z, y, x].
    #     """
    #     pass

    # @abstractmethod
    # def get_grid_x(self, long_var_name):
    #     """
    #     Only return something for variables with a rectilinear, structured or unstructured grid. Otherwise raise ValueError.
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     Numpy array of doubles: x coordinate of grid cell center for each grid cell, in the same order as the values returned by function get_value.
    #                      For a rectilinear grid: x coordinate of column center for each column.
    #     """
    #     pass

    # @abstractmethod
    # def get_grid_y(self, long_var_name):
    #     """
    #     Only return something for variables with a rectilinear, structured or unstructured grid. Otherwise raise ValueError.
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     Numpy array of doubles: y coordinate of grid cell center for each grid cell, in the same order as the values returned by function get_value.
    #                      For a rectilinear grid: y coordinate of row center for each row.
    #     """
    #     pass

    # @abstractmethod
    # def get_grid_z(self, long_var_name):
    #     """
    #     Only return something for variables with a rectilinear, structured or unstructured grid. Otherwise raise ValueError.
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     Numpy array of doubles: z coordinate of grid cell center for each grid cell, in the same order as the values returned by function get_value.
    #                      For a rectilinear grid: z coordinate of layer center for each layer.
    #     """
    #     pass

    # @abstractmethod
    # def get_grid_connectivity(self, long_var_name):
    #     """
    #     Only return something for variables with an unstructured grid. Otherwise raise ValueError.
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     ???
    #     """
    #     pass

    # @abstractmethod
    # def get_grid_offset(self, long_var_name):
    #     """
    #     Only return something for variables with an unstructured grid. Otherwise raise ValueError.
    #     Input parameters:
    #     String long_var_name: identifier of a variable in the model.
    #     Return value:
    #     ???
    #     """
    #     pass
    
    
class EBmi(Bmi):
    
    @abstractmethod
    def initialize_config(self, config_fn):
        """
        First step of two-phase initialize. In this step only the configuration is read in.
        This allows a user to then change settings and parameters before fully initializing the model
        """
        pass
    
    @abstractmethod
    def initialize_model(self, source_directory):
        """
        Second step of two-phase initialize. In this step the model is setup, and will now allow
        reading/setting values.
        """
        pass
    
    @abstractmethod
    def set_start_time(self, start_time):
        """
        Set the start time of the model. Can usually only be called after initialize_config
        and before initialize_model. Expects a value in the time units of the model
        """
        pass
    
    @abstractmethod
    def set_end_time(self, end_time):
        """
        Set the end time of the model. Can usually only be called after initialize_config
        and before initialize_model. Expects a value in the time units of the model.
        """
        pass
    
    @abstractmethod
    def get_attribute_names(self):
        """
        Gets a list of all supported attributes for this model. Attributes can be considered
        the meta-data of a model, for instance author, version, model specific settings, etc.
        """
        pass
    
    @abstractmethod
    def get_attribute_value(self, attribute_name):
        """
        Gets the value of a certain attribute for this model. Attributes can be considered
        the meta-data of a model, for instance author, version, model specific settings, etc.
        """
        pass
    
    @abstractmethod
    def set_attribute_value(self, attribute_name, attribute_value):
        """
        Sets the value of a certain attribute for this model. Usually only string values are allowed.
        Attributes can be considered the meta-data of a model, for instance author, version, model specific settings, etc.
        """
        pass
    
    # @abstractmethod
    # def save_state(self, destination_directory):
    #     """
    #     Ask the model to write its complete internal current state to one or more state files in the given directory.
    #     Afterwards the given directory should only contain the state files and nothing else.
    #     Input parameters:
    #     File destination_directory: the directory in which the state files should be written.
    #     """
    #     pass
    
    # @abstractmethod
    # def load_state(self, source_directory):
    #     """
    #     Ask the model to load its complete internal current state from one or more state files in the given directory.
 
    #     Input parameters:
    #     File source_directory: the directory from which the state files should be read.
    #     """
    #     pass



class GBmi(EBmi):
    """
    Variable Information Functions
    """

    def get_var_shape(self, long_var_name):
        """
        Input parameters:
        String long_var_name: identifier of a variable in the model.
        Return value:
        Tuple of Integer: Variable shape
        """
        return self.get_value(long_var_name).shape

    """
    Grid Information Functions
    """

    def get_grid_type(self):
        """
        Return value:
        GridType type of bmi model 
        """
        return self.get_grid().get_type()

    @abstractmethod
    def get_grid(self):
        """
        Create Grid definition object, see grid.py

        Return value:
        grid.py object
        """
        pass


    """
    Model control function
    """
    @abstractmethod
    def spinup(self):
        """PCR specific spinup function"""
        pass


    """
    set and get attribute / config 
    """
    @abstractmethod
    def write_config(self):
        """write adapted config to file. just before initializing
        only for models which do not allow for direct access to model config via bmi"""
        pass

    @abstractmethod
    def set_out_dir(self, out_dir):
        pass
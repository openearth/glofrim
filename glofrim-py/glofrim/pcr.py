import numpy as np
import sys
from configparser import ConfigParser
import logging
from os.path import join, isfile, abspath
from datetime import datetime, date
import rasterio

from utils import setlogger
from gbmi import GBmi, BmiGridType
from gbmi_lib import configread

class PCR(GBmi):
    """
    Glofrim implementation of the PCR BMI adaptor.
    """
    _name = 'PCRGLOB-WB'
    _var_units = {'var1': 'unit1'}
    _input_var_names = ['var1']
    _output_var_names = ['var1']

    def __init__ (self):
        # import original PCR bmi 
        from pcrglobwb_bmi_v203 import pcrglobwb_bmi as _bmi
        self._bmi = _bmi.pcrglobwbBMI()
        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False

    """
    Model Control Functions
    """
    def initialize_config(self, config_fn):
        cf = ConfigParser(inline_comment_prefixes=('#'))
        self._config_fn = config_fn
        self._config = configread(config_fn, encoding='utf-8', cf=cf)
        # model time
        self._dt = 1. 
        self._timeunit = 'day'
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model files
        self._ddir = abspath(self.get_attribute_value('globalOptions:inputDir'))
        self._lm_fn = join(self._ddir, self.get_attribute_value('globalOptions:landmask'))
        self._ldd_fn = join(self._ddir, self.get_attribute_value('routingOptions:lddMap'))
        if not isfile(self._lm_fn): raise IOError('landmask file not found')
        if not isfile(self._ldd_fn): raise IOError('ldd file not found')

    def initialize_model(self, source_directory=None):
        #NOTE: PCR does not use a source_directory
        if not hasattr(self, '_config_fn'):
            raise Warning('run initialize_config before initialize_model')
        self._bmi.initialize(self._config_fn)
        self.initialized = True

    def initialize(self, config_fn):
        self.initialize_config(config_fn)
        self.initialize_model()
        self.logger.info('Model initialized')
    
    def update(self):
        # if self._t >= self._endTime:
		#     raise Exception("endTime already reached, model not updated")
        self._bmi.update(dt=self._dt)
        self._t = self.get_current_time()

    def update_until (self, t):
        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update ()

    def finalize(self):
        self._bmi.finalize()


    """
    Model Information Functions
    """

    def get_component_name (self):
        return self._name

    def get_input_var_names (self):
        return self._input_var_names

    def get_output_var_names (self):
        return self._output_var_names


    """
    Variable Information Functions
    """

    def get_var_type (self, long_var_name):
        return str(self.get_value(long_var_name).dtype)

    def get_var_units (self, long_var_name):
        return self._var_units[long_var_name]

    def get_var_rank (self, long_var_name):
        return self.get_value(long_var_name).ndim

    def get_var_size(self, long_var_name):
        return self.get_value(long_var_name).size

    def get_var_nbytes(self, long_var_name):
        return self.get_value(long_var_name).nbytes

    
    def get_start_time(self):
        if self.initialized:
            startTime = self._bmi.get_start_time()
        else:
            startTime = self.get_attribute_value('globalOptions:startTime')
            startTime = datetime.strptime(startTime, '%Y-%m-%d').date()
        self._startTime = startTime
        return self._startTime
    
    def get_current_time(self):
        return self._bmi.get_current_time()

    def get_end_time(self):
        if self.initialized:
            endTime = self._bmi.get_end_time()
        else:
            endTime = self.get_attribute_value('globalOptions:endTime')
            endTime = datetime.strptime(endTime, '%Y-%m-%d').date()
        self._endTime = endTime
        return self._endTime

    def get_time_step(self):
        #NOTE get_time_step in PCR bmi returns timestep as int instead of dt
        return self._dt 
        
    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, fill_value=-999):
        # additional fill_value argument required to translate pcr maps to numpy arrays
        return np.asarray(self._bmi.get_var(long_var_name, missingValues=fill_value))

    def get_value_at_indices(self, long_var_name, inds, fill_value=-999):
        # always use 1d inds
        if not ((isinstance(inds, np.ndarray)) and (inds.ndim == 1)):
            raise ValueError('indices should be 1d arrays')
        return self.get_value(long_var_name, fill_value=fill_value).flat[inds]

    def set_value(self, long_var_name, src, fill_value=-999):
        self._bmi.set_var(long_var_name, src, missingValues=fill_value)

    def set_value_at_indices(self, long_var_name, inds, src, fill_value=-999):
        # always use 1d inds
        if not ((isinstance(inds, np.ndarray)) and (inds.ndim == 1)):
            raise ValueError('indices should be 1d arrays')
        val = self.get_value(long_var_name, missingValues=fill_value)
        val.flat[inds] = src
        self._bmi.set_var(long_var_name, val, missingValues=fill_value)
        #TODO do we still need to set data or is get_value a pointer?

    
    """
    Grid Information Functions
    """
    
    def get_grid_type(self, long_var_name):
        return BmiGridType.RECTILINEAR

    def get_grid_transform(self):
        if not hasattr(self, 'grid_transform'):
            self.logger.info('Getting PCR model transform based on landmask file')
            with rasterio.open(self._lm_fn, 'r') as ds:
                self.grid_transform = ds.transform
        return self.grid_transform

    def get_grid_shape(self, **kwargs):
        if not hasattr(self, 'grid_shape'):
            self.logger.info('Getting PCR model transform based on landmask file')
            with rasterio.open(self._lm_fn, 'r') as ds:
                self.grid_shape = ds.shape
        return self.grid_shape

    def grid_index(self, x, y, **kwargs):
        """Get PCR indices (1d) of xy coordinates."""
        if not hasattr(self, '_grid_index'):
            self.logger.info('Getting PCR model index based on landmask file')
            with rasterio.open(self._lm_fn, 'r') as ds:
                self._grid_index = ds.index
        if not hasattr(self, '_mask'):
            self.logger.info('Getting PCR model mask based on landmask file')
            with rasterio.open(self._lm_fn, 'r') as ds:
                self._mask = ds.read(1)==1
        # get row, cols
        self.logger.info('Getting PCR model indices of xy coordinates.')
        r, c = self._grid_index(x, y)
        r, c = np.asarray(r).astype(int), np.asarray(c).astype(int)
        # check if inside domain
        nrows, ncols = self.get_grid_shape()
        inside = np.logical_and.reduce((r>=0, r<nrows, c>=0, c<ncols))
        # get valid cells (inside landmask and mask)
        valid = np.zeros_like(inside, dtype=bool) 
        inds = np.ones_like(inside, dtype=int)*-1 
        # import pdb; pdb.set_trace()
        if np.any(inside):
            valid[inside] = self._mask[r[inside], c[inside]]
        # calculate 1d index 
        # NOTE invalid indices have value -1
        if np.any(valid):
            inds[valid] = np.ravel_multi_index((r[valid], c[valid]), (nrows, ncols))
        return inds

    def get_drainage_direction(self):
        from nb.nb_io import read_dd_pcraster
        # read file with pcr readmap
        if not hasattr(self, 'grid_transform'):
            self.get_grid_transform()
        self.dd = read_dd_pcraster(self._ldd_fn, self.grid_transform, nodata=255)


    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, (datetime, date)):
            start_time = start_time.strftime("%Y-%m-%d")
        try:
            start_time.strftime("%Y-%m-%d") 
        except ValueError:
            raise ValueError('wrong date format, use "yyyy-mm-dd"')

        self.set_attribute_value('globalOptions:startTime', start_time)
       
    def set_end_time(self, end_time):
        if isinstance(end_time, (datetime, date)):
            end_time = end_time.strftime("%Y-%m-%d")
        try:
            end_time.strftime("%Y-%m-%d") 
        except ValueError:
            raise ValueError('wrong date format, use "yyyy-mm-dd"')
        self.set_attribute_value('globalOptions:endTime', end_time)

    def get_attribute_names(self):
        attr = []
        for sect in self._config.sections():
            for opt in self._config.options(sect):
                attr.append(sect + ":" + opt)
        return attr
    
    def get_attribute_value(self, attribute_name):
        self.logger.debug("get_attribute_value: " + attribute_name)
        attrpath = attribute_name.split(":")
        if len(attrpath) == 2:
            return self._config.get(attrpath[0], attrpath[1])
        else:
            msg = "Attributes should follow the name:option convention"
            self.logger.warn(msg)
            raise Warning(msg)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        self.logger.debug("set_attribute_value: " + attribute_value)
        attrpath = attribute_name.split(":")
        if len(attrpath) == 2:
            self._config.set(attrpath[0], attrpath[1], attribute_value)
        else:
            msg = "Attributes should follow the name:option convention"
            self.logger.warn(msg)
            raise Warning(msg)


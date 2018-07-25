import numpy as np
import sys
from configparser import ConfigParser
import logging
import os
from os.path import join, isfile, abspath, dirname, basename, normpath
from datetime import datetime, timedelta
import rasterio

from utils import setlogger
from gbmi import GBmi
from grids import RGrid
import glofrim_lib as glib 


class PCR(GBmi):
    """
    csdms BMI implementation of the PCR BMI adaptor for GLOFRIM.
    """
    _name = 'PCR'
    _long_name = 'PCRGLOB-WB'
    _version = '2.0.3'
    _var_units = {'runoff': 'm/day', 'discharge': 'm3/s', 'cellArea': 'm2'}
    _input_var_names = ['cellArea']
    _output_var_names = ['runoff', 'discharge']
    _area_var_name = 'cellArea'
    _timeunit = 'days'
    _dt = timedelta(days=1) # NOTE: this is not an optoins in PCR

    def __init__(self):
        # import original PCR bmi 
        from pcrglobwb_bmi_v203 import pcrglobwb_bmi as _bmi
        self._bmi = _bmi.pcrglobwbBMI()
        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False
        self.grid = None

    """
    Model Control Functions
    """
    def initialize_config(self, config_fn):
        # config settings
        if self.initialized:
            raise Warning("model already initialized, it's therefore not longer possible to initialize the config")
        self._config_fn = abspath(config_fn)
        self._config = glib.configread(self._config_fn, encoding='utf-8', 
                                       cf=ConfigParser(inline_comment_prefixes=('#')))
        self._datefmt = "%Y-%m-%d"
        # model time
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model files
        self._outdir= abspath(self.get_attribute_value('globalOptions:outputDir'))
        self.logger.info('Config initialized')

    def initialize_model(self, **kwargs):
        if not hasattr(self, '_config_fn'):
            raise Warning('run initialize_config before initialize_model')
        self.write_config() # write updated config to file as bmi does not allow direct access
        self._bmi.initialize(self._config_fn)
        self.initialized = True
        self.logger.info('Model initialized')
        # reset model time to make sure it is consistent with the model
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime

    def initialize(self, config_fn):
        self.initialize_config(config_fn)
        self.initialize_model()
            
    def update(self):
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        self._bmi.update(dt=self._dt.days)
        self._t += self._dt
        self.logger.info('updated model to datetime {}'.format(self._t.strftime("%Y-%m-%d")))

    def update_until(self, t):
        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update()

    def spinup(self):
        """PCR specific spinup function"""
        self._bmi.spinup()

    def finalize(self):
        self._bmi.finalize()


    """
    Variable Information Functions
    """
    
    def get_start_time(self):
        if self.initialized:
            # date to datetime object
            startTime = datetime.combine(self._bmi.get_start_time(), datetime.min.time())
        else:
            startTime = self.get_attribute_value('globalOptions:startTime')
            startTime = datetime.strptime(startTime, self._datefmt)
        self._startTime = startTime
        return self._startTime
    
    def get_current_time(self):
        if self.initialized:
            return self._t
        else:
            return self.get_start_time()

    def get_end_time(self):
        if self.initialized:
            # date to datetime object
            endTime = datetime.combine(self._bmi.get_end_time(), datetime.min.time())
        else:
            endTime = self.get_attribute_value('globalOptions:endTime')
            endTime = datetime.strptime(endTime, self._datefmt)
        self._endTime = endTime
        return self._endTime

    def get_time_step(self):
        return self._dt 
        
    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, fill_value=-999, **kwargs):
        # additional fill_value argument required to translate pcr maps to numpy arrays
        var = np.asarray(self._bmi.get_var(long_var_name, missingValues=fill_value, **kwargs))
        var = np.where(var == fill_value, np.nan, var)
        return var

    def get_value_at_indices(self, long_var_name, inds, fill_value=-999, **kwargs):
        return self.get_value(long_var_name, fill_value=fill_value, **kwargs).flat[inds]

    def set_value(self, long_var_name, src, fill_value=-999, **kwargs):
        src = np.where(np.isnan(src), fill_value, src).astype(self.get_var_type(long_var_name))
        self._bmi.set_var(long_var_name, src, missingValues=fill_value, **kwargs)

    def set_value_at_indices(self, long_var_name, inds, src, fill_value=-999,  **kwargs):
        val = self.get_value(long_var_name, fill_value=fill_value, **kwargs)
        val.flat[inds] = src
        self.set_value(long_var_name, val, **kwargs)
        

    """
    Grid Information Functions
    """
    def get_grid(self):

        if not hasattr(self, 'grid') or (self.grid is None):
            # ldd is used for coupling to routing / flood model
            _indir = abspath(self.get_attribute_value('globalOptions:inputDir'))
            _ldd_fn = glib.getabspath(self.get_attribute_value('routingOptions:lddMap'), _indir)
            if not isfile(_ldd_fn): raise IOError('ldd file not found')
            # landmask used for masking out coordinates outside mask
            _lm_fn = glib.getabspath(self.get_attribute_value('globalOptions:landmask'), _indir)
            if not isfile(_lm_fn): raise IOError('landmask file not found')
            self.logger.info('Getting rgrid info based on {}'.format(basename(_lm_fn)))
            with rasterio.open(_lm_fn, 'r') as ds:
                self.grid = RGrid(ds.transform, ds.height, ds.width, crs=ds.crs, mask=ds.read(1)==1)

            # read file with pcr readmap
            self.logger.info('Getting drainage direction from {}'.format(basename(_ldd_fn)))
            self.grid.set_dd(_ldd_fn, ddtype='ldd', nodata=255)
        return self.grid

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, datetime):
            start_time = start_time.strftime(self._datefmt)
        try:
            self._startTime = datetime.strptime(start_time, self._datefmt) 
            self._t = self._startTime
        except ValueError:
            raise ValueError('wrong date format, use "yyyy-mm-dd"')
        self.set_attribute_value('globalOptions:startTime', start_time)
       
    def set_end_time(self, end_time):
        if isinstance(end_time, datetime):
            end_time = end_time.strftime(self._datefmt)
        try:
            self._endTime = datetime.strptime(end_time, self._datefmt) 
        except ValueError:
            raise ValueError('wrong date format, use "yyyy-mm-dd"')
        self.set_attribute_value('globalOptions:endTime', end_time)

    def set_out_dir(self, out_dir):
        self.set_attribute_value('globalOptions:outputDir', abspath(out_dir))

    def get_attribute_names(self):
        glib.configcheck(self, self.logger)
        return glib.configattr(self._config)
    
    def get_attribute_value(self, attribute_name):
        glib.configcheck(self, self.logger)
        self.logger.debug("get_attribute_value: " + attribute_name)
        return glib.configget(self._config, attribute_name)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        glib.configcheck(self, self.logger)
        self.logger.debug("set_attribute_value: " + attribute_value)
        return glib.configset(self._config, attribute_name, attribute_value)

    def write_config(self):
        """write adapted config to file. just before initializing
        only for models which do not allow for direct access to model config via bmi"""
        self._config_fn = glib.write_config(self, self._config, self._config_fn, self.logger)
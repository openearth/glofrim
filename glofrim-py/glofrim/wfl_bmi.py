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


class WFL(GBmi):
    """
    csdms BMI implementation of the WFLOW BMI adaptor for GLOFRIM.
    """
    _name = 'WFL'
    _long_name = 'wflow'
    _version = ''
    _var_units = {}
    _input_var_names = []
    _output_var_names = []
    _area_var_name = 'csize'
    _timeunit = 'seconds'

    def __init__(self):
        # import original PCR bmi 
        import wflow.wflow_bmi as _bmi
        self._bmi = _bmi.wflowbmi_csdms()
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
        self.logger.info('Read ini at {}'.format(self._config_fn))
        self._bmi.initialize_config(str(self._config_fn))
        self._config = self._bmi.config
        self._datefmt = "%Y-%m-%d %H:%M:%S"
        if self.get_attribute_value('run:runlengthdetermination') == 'steps':
            raise Warning('set WFLOW run:runlengthdetermination to "intervals" in order to use the same time convention as in GLOFRIM')
        # model time
        self._dt = self.get_time_step()
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model files
        self.logger.info('Config initialized')

        namesroles = self._bmi.dynModel.wf_supplyVariableNamesAndRoles()
        for n,r,u in namesroles:
            if r == 0 or r == 3: #input  or parameter
                self._input_var_names.append(n)
            elif r == 1 or r == 2: # output or state
                self._output_var_names.append(n)
            self._var_units[n] = u

    def initialize_model(self, **kwargs):
        if not hasattr(self, '_config_fn'):
            raise Warning('run initialize_config before initialize_model')
        self._bmi.initialize_model()
        self.initialized = True
        self.logger.info('Model initialized')
        # reset model time to make sure it is consistent with the model
        self._dt = self.get_time_step()
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime

    def initialize(self, config_fn):
        self.initialize_config(config_fn)
        self.initialize_model()
            
    def update(self):
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        self._bmi.update()
        self._t += self.get_time_step()

    def update_until(self, t):
        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update()

    def spinup(self):
        """PCR specific spinup function"""
        raise NotImplementedError()

    def finalize(self):
        self._bmi.finalize()


    """
    Variable Information Functions
    """
    
    def get_start_time(self):
        # if self.initialized:
        #     startTime = datetime.utcfromtimestamp(self._bmi.get_start_time())
        # else:
        #     startTime = self.get_attribute_value('run:starttime')
        #     startTime = datetime.strptime(startTime, self._datefmt)
        startTime = datetime.utcfromtimestamp(self._bmi.get_start_time())
        self._startTime = startTime
        return self._startTime
    
    def get_current_time(self):
        if self.initialized:
            #
            return datetime.utcfromtimestamp(self._bmi.get_current_time())
        else:
            return self.get_start_time()

    def get_end_time(self):
        # if self.initialized:
        #     endTime = datetime.utcfromtimestamp(self._bmi.get_end_time())
        # else:
        #     endTime = self.get_attribute_value('run:endtime')
        #     endTime = datetime.strptime(endTime, self._datefmt)
        endTime = datetime.utcfromtimestamp(self._bmi.get_end_time())        
        self._endTime = endTime
        return self._endTime

    def get_time_step(self):
        if not hasattr(self, '_dt'):
            self._dt = timedelta(**{self.get_time_units(): self._bmi.get_time_step()})
        return self._dt 
        
    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, **kwargs):
        # additional fill_value argument required to translate pcr maps to numpy arrays
        return np.asarray(self._bmi.get_value(long_var_name))

    def get_value_at_indices(self, long_var_name, inds, **kwargs):
        return self.get_value(long_var_name).flat[inds]

    def set_value(self, long_var_name, src, **kwargs):
        self._bmi.set_value(long_var_name, src.astype(self.get_var_type(long_var_name)))

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        val = self.get_value(long_var_name)
        val.flat[inds] = src
        self.set_value(long_var_name, val)
        

    """
    Grid Information Functions
    """
    def get_grid(self, mapdir=None):
        if not hasattr(self, 'grid') or (self.grid is None):
            if mapdir is None:
                mapdir = join(dirname(self._config_fn), 'staticmaps')
            # ldd is used for coupling to routing / flood model
            _ldd_fn = join(mapdir, 'wflow_ldd.map')
            if not isfile(_ldd_fn): raise IOError('ldd file not found')
            # landmask used for masking out coordinates outside mask
            _lm_fn = join(mapdir, 'wflow_subcatch.map')
            if not isfile(_lm_fn): raise IOError('subcatch file not found')
            self.logger.info('Getting rgrid info based on {}'.format(basename(_lm_fn)))
            with rasterio.open(_lm_fn, 'r') as ds:
                self.grid = RGrid(ds.transform, ds.height, ds.width, crs=ds.crs, mask=ds.read(1)>=1)
            # read file with pcr readmap
            self.logger.info('Getting drainage direction from {}'.format(basename(_ldd_fn)))
            self.grid.set_dd(_ldd_fn, ddtype='ldd', nodata=255)
        return self.grid

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, str):
            start_time = datetime.strptime(start_time, self._datefmt) 
        t = (start_time - datetime.utcfromtimestamp(0))
        self._startTime = start_time
        self._t = start_time
        self._bmi.set_start_time(t.days*86400+t.seconds)
       
    def set_end_time(self, end_time):
        if isinstance(end_time, str):
            end_time = datetime.strptime(end_time, self._datefmt) 
        t = (end_time - datetime.utcfromtimestamp(0))
        self._endTime = end_time
        self._bmi.set_end_time(t.days*86400+t.seconds)

    def set_out_dir(self, out_dir):
        # NOTE this needs te be set somewhere pre initialize config
        # the file setup is hardcoded in wflow
        pass
    
    def get_attribute_names(self):
        glib.configcheck(self, self.logger)
        return glib.configattr(self._config)
    
    def get_attribute_value(self, attribute_name):
        glib.configcheck(self, self.logger)
        self.logger.debug("get_attribute_value: " + attribute_name)
        return self._bmi.get_attribute_value(attribute_name)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        glib.configcheck(self, self.logger)
        self.logger.info("set_attribute_value {} -> {} ".format(attribute_name, attribute_value))
        # update config
        glib.configset(self._config, attribute_name, attribute_value)
        # set wfl config
        return self._bmi.set_attribute_value(attribute_name, attribute_value)

    def write_config(self):
        """write adapted config to file. just before initializing
        only for models which do not allow for direct access to model config via bmi"""
        glib.write_config(self, self._config, self._config_fn, self.logger)
import numpy as np
import sys
from configparser import ConfigParser
import logging
import os
from os.path import join, isfile, abspath, dirname, basename, normpath
from datetime import datetime, timedelta
import rasterio

from utils import setlogger
from gbmi import GBmi, GBmiGridType, GBmiModelType
import gbmi_lib as glib 


class PCR(GBmi):
    """
    csdms BMI implementation of the PCR BMI adaptor for GLOFRIM.
    """
    _name = 'PCRGLOB-WB'
    _version = '2.0.3'
    _var_units = {'runoff': 'm/day', 'discharge': 'm3/s', 'cellArea': 'm2'}
    _input_var_names = ['var1']
    _output_var_names = ['runoff', 'discharge']
    _area_var_name = 'cellArea'

    def __init__(self):
        # import original PCR bmi 
        from pcrglobwb_bmi_v203 import pcrglobwb_bmi as _bmi
        self._bmi = _bmi.pcrglobwbBMI()
        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False

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
        self._dt = timedelta(days=1)
        self._timeunit = 'day'
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model files
        self._indir = abspath(self.get_attribute_value('globalOptions:inputDir'))
        self._outdir= abspath(self.get_attribute_value('globalOptions:outputDir'))
        # ldd is used for coupling to routing / flood model
        self._ldd_fn = glib.getabspath(self.get_attribute_value('routingOptions:lddMap'), self._indir)
        if not isfile(self._ldd_fn): raise IOError('ldd file not found')
        # landmask used for masking out coordinates outside mask
        self._lm_fn = glib.getabspath(self.get_attribute_value('globalOptions:landmask'), self._indir)
        if not isfile(self._lm_fn): raise IOError('landmask file not found')
        self.logger.info('Config initialized')

    def initialize_model(self, **kwargs):
        if not hasattr(self, '_config_fn'):
            raise Warning('run initialize_config before initialize_model')
        self.write_config() # write updated config to file as bmi does not allow direct access
        self._bmi.initialize(self._config_fn)
        self.initialized = True
        self.logger.info('Model initialized')

    def initialize(self, config_fn):
        self.initialize_config(config_fn)
        self.initialize_model()
            
    def update(self):
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        self._bmi.update(dt=self.get_time_step().days)
        self._t += self.get_time_step()

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
    Model Information Functions
    """
    def get_model_type(self):
        return GBmiModelType.HYDROLOGICAL 

    def get_component_name(self):
        return self._name

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names


    """
    Variable Information Functions
    """

    def get_var_type(self, long_var_name):
        return str(self.get_value(long_var_name).dtype)

    def get_var_units(self, long_var_name):
        return self._var_units[long_var_name]

    def get_var_rank(self, long_var_name):
        return self.get_value(long_var_name).ndim

    def get_var_size(self, long_var_name):
        return self.get_value(long_var_name).size

    def get_var_shape(self, long_var_name):
        return self.get_value(long_var_name).shape

    def get_var_nbytes(self, long_var_name):
        return self.get_value(long_var_name).nbytes

    
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

    def get_drainage_direction(self):
        from nb.nb_io import read_dd_pcraster
        # read file with pcr readmap
        if not hasattr(self, 'grid_transform'):
            self.get_grid_transform()
        self.logger.info('Getting drainage direction from {}'.format(basename(self._ldd_fn)))
        self.dd = read_dd_pcraster(self._ldd_fn, self.grid_transform, nodata=255)

    """
    Grid Information Functions
    """
    
    def get_grid_type(self):
        return GBmiGridType.RECTILINEAR

    def get_grid_transform(self):
        if not hasattr(self, 'grid_transform'):
            self.logger.info('Getting grid transform based on {}'.format(basename(self._lm_fn)))
            with rasterio.open(self._lm_fn, 'r') as ds:
                self.grid_transform = ds.transform
        return self.grid_transform

    def get_grid_shape(self):
        if not hasattr(self, 'grid_shape'):
            self.logger.info('Getting grid shape based on {}'.format(basename(self._lm_fn)))
            with rasterio.open(self._lm_fn, 'r') as ds:
                self.grid_shape = ds.shape
        return self.grid_shape

    def get_grid_bounds(self):
        if not hasattr(self, 'grid_bounds'):
            self.logger.info('Getting grid bounds based on {}'.format(basename(self._lm_fn)))
            with rasterio.open(self._lm_fn, 'r') as ds:
                self.grid_bounds = ds.bounds
        return self.grid_bounds

    def get_grid_res(self):
        if not hasattr(self, 'grid_res'):
            self.logger.info('Getting grid bounds based on {}'.format(basename(self._lm_fn)))
            with rasterio.open(self._lm_fn, 'r') as ds:
                self.grid_res = ds.res
        return self.grid_res

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

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, datetime):
            start_time = start_time.strftime(self._datefmt)
        try:
            datetime.strptime(start_time, self._datefmt) 
        except ValueError:
            raise ValueError('wrong date format, use "yyyy-mm-dd"')

        self.set_attribute_value('globalOptions:startTime', start_time)
       
    def set_end_time(self, end_time):
        if isinstance(end_time, datetime):
            end_time = end_time.strftime(self._datefmt)
        try:
            datetime.strftime(end_time, self._datefmt) 
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
        glib.configcheck(self, self.logger)
        out_dir = dirname(self._config_fn)
        bname = basename(self._config_fn).split('.')
        new_fn = '{}_glofrim.{}'.format('.'.join(bname[:-1]), bname[-1])
        self._config_fn = join(out_dir, new_fn)
        if isfile(self._config_fn):
            os.unlink(self._config_fn)
            self.logger.warn("{:s} file overwritten".format(self._config_fn))
        glib.configwrite(self._config, self._config_fn, encoding='utf-8')
        self.logger.info('Ini file written to {:s}'.format(self._config_fn))
import numpy as np
import sys
from configparser import ConfigParser
import logging
import os
from os.path import join, isfile, abspath, dirname, basename, relpath
from datetime import datetime, timedelta
import rasterio
import re

from bmi.wrapper import BMIWrapper as _bmi

from utils import setlogger
from gbmi import GBmi 
from grids import UCGrid
import glofrim_lib as glib 

class CMF(GBmi):
    """
    Glofrim implementation of the PCR BMI adaptor.
    """
    _name = 'CMF'
    _long_name = 'CaMa-Flood'
    _version = '3.6.2'
    _var_units = {'roffin': 'm', # NOTE: unit set by user, this is edited when initialize config
                  'outflw': 'm3/s', 'runoff': 'm3/s'}
    _input_var_names = ['roffin']
    _output_var_names = ['outflw']
    _area_var_name = ''
    _timeunit = 'seconds'

    def __init__(self, engine):
        self._bmi = _bmi(engine = engine)
        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False
        self.grid = None

    """
    Model Control Functions
    """
    def initialize_config(self, config_fn):
        if self.initialized:
            raise Warning("model already initialized, it's therefore not longer possible to initialize the config")
        # config settings
        if basename(config_fn) == 'input_flood.nam':
            msg = '"input_flood.nam" is protected and not allowed to be used as configuration file name for \
                    CaMa-Flood. Rename the input configuration file to e.g. "input_flood.temp"'
            raise Warning(msg)
        self._config_fn = abspath(config_fn)
        self._config = glib.configread(self._config_fn, encoding='utf-8', cf=NamConfigParser())
        # model time
        self._dt = self.get_time_step()
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model units
        rofunit = float(self.get_attribute_value('CONF:DROFUNIT').replace('D', 'e'))
        self._var_units['roffin'] = {1: 'm', 1e-3:'mm'}[rofunit]
        # model files
        _root = dirname(self._config_fn)
        # mapdir where nextxy data is found
        self._mapdir = dirname(glib.getabspath(str(self.get_attribute_value('MAP:CNEXTXY').strip('"')), _root))
        self._outdir = glib.getabspath(str(self.get_attribute_value('OUTPUT:COUTDIR').strip('"')), _root)
        self.logger.info('Config initialized')

    def initialize_model(self, **kwargs):
        if not hasattr(self, '_config_fn'):
            raise Warning('Run initialize_config before initialize_model')
        self.write_config() # write updated config to file as bmi does not allow direct access
        self._bmi.initialize(self._config_fn)
        self.initialized = True
        self.logger.info('Model initialized')

    def initialize(self, config_fn):
        if not hasattr(self, '_config'):
            self.initialize_config(config_fn)
        self.initialize_model()
            
    def update(self):
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        self._bmi.update(dt=self.get_time_step().days)
        self._t += self.get_time_step()
        self.logger.info('updated model to datetime {}'.format(self._t.strftime("%Y-%m-%d")))

    def update_until(self, t):
        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update()

    # not defined in CMF
    def spinup(self):
        """PCR specific spinup function"""
        raise NotImplementedError()

    def finalize(self):
        self._bmi.finalize()


    """
    Variable Information Functions
    """
    
    def get_start_time(self):
        if self.initialized:
            startTime = CMFtime_2_datetime(self._bmi.get_start_time())
        else:
            yr = int(self.get_attribute_value('SIMTIME:ISYEAR'))
            m = int(self.get_attribute_value('SIMTIME:ISMON'))
            d = int(self.get_attribute_value('SIMTIME:ISDAY'))
            startTime = datetime(yr, m, d, 0, 0)
        self._startTime = startTime
        return self._startTime
    
    def get_current_time(self):
        if self.initialized:
            return CMFtime_2_datetime(self._bmi.get_current_time())
        else:
            return self.get_start_time()

    def get_end_time(self):
        if self.initialized:
            endTime = CMFtime_2_datetime(self._bmi.get_end_time())
        else:
            yr = int(self.get_attribute_value('SIMTIME:IEYEAR'))
            m = int(self.get_attribute_value('SIMTIME:IEMON'))
            d = int(self.get_attribute_value('SIMTIME:IEDAY'))
            endTime = datetime(yr, m, d, 0, 0)
        self._endTime = endTime
        return self._endTime

    def get_time_step(self):
        if not hasattr(self, '_dt'):
            dt = int(self.get_attribute_value('CONF:DT'))
            self._dt = timedelta(**{self.get_time_units(): dt})
        return self._dt 
        
    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, **kwargs):
        return np.asarray(self._bmi.get_var(long_var_name))

    def get_value_at_indices(self, long_var_name, inds, **kwargs):
        return self.get_value(long_var_name).flat[inds]

    def set_value(self, long_var_name, src, **kwargs):
        self._bmi.set_var(long_var_name, src)

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        val = self.get_value(long_var_name)
        val.flat[inds] = src
        self._bmi.set_var(long_var_name, val)

    """
    Grid Information Functions
    """

    def get_grid(self, fn_catmxy=r'./hires/reg.catmxy.tif', **kwargs):
        if not hasattr(self, 'grid') or (self.grid is None):
            _root = dirname(self._config_fn)
            _nextxy_fn = glib.getabspath(str(self.get_attribute_value('MAP:CNEXTXY').strip('"')), _root)
            _nextxy_fn = _nextxy_fn.replace('.bin', '.tif')
            if not isfile(_nextxy_fn): raise IOError('nextxy file {} not found'.format(_nextxy_fn))
            # fn_catmx fn may be relative to folder where nextxy is found
            fn_catmxy = glib.getabspath(fn_catmxy, dirname(_nextxy_fn))
            if not isfile(fn_catmxy):
                raise IOError("{} file not found".format(fn_catmxy))
            self.logger.info('Getting Unit-Catchment Grid info based on {}'.format(basename(_nextxy_fn)))
            with rasterio.open(_nextxy_fn, 'r') as ds:
                self.grid = UCGrid(ds.transform, ds.height, ds.width, fn_catmxy=fn_catmxy, crs=ds.crs)
            # set drainage direction
            self.logger.info('Getting drainage direction from {}'.format(basename(_nextxy_fn)))
            self.grid.set_dd(_nextxy_fn, ddtype='nextxy', **kwargs)
        return self.grid

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, str):
            start_time = datetime.strptime(start_time, "%Y-%m-%d") 
        if not isinstance(start_time, datetime):
            raise ValueError("invalid date type")
        self.set_attribute_value('SIMTIME:ISYEAR', start_time.year)
        self.set_attribute_value('SIMTIME:ISMON', start_time.month)
        self.set_attribute_value('SIMTIME:ISDAY', start_time.day)
       
    def set_end_time(self, end_time):
        if isinstance(end_time, str):
            end_time = datetime.strptime(end_time, "%Y-%m-%d") 
        if not isinstance(end_time, datetime):
            raise ValueError("invalid date type")
        self.set_attribute_value('SIMTIME:IEYEAR', end_time.year)
        self.set_attribute_value('SIMTIME:IEMON', end_time.month)
        self.set_attribute_value('SIMTIME:IEDAY', end_time.day)

    def set_out_dir(self, out_dir):
        self.set_attribute_value('OUTPUT:COUTDIR', abspath(out_dir))

    def get_attribute_names(self):
        glib.configcheck(self, self.logger)
        return glib.configattr(self._config)
    
    def get_attribute_value(self, attribute_name):
        glib.configcheck(self, self.logger)
        self.logger.debug("get_attribute_value: {}".format(attribute_name))
        return glib.configget(self._config, attribute_name)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        glib.configcheck(self, self.logger)
        self.logger.debug("set_attribute_value: {} -> {}".format(attribute_name, attribute_value))
        return glib.configset(self._config, attribute_name, str(attribute_value))

    def set_inpmat_attrs(self):
        # set new inpmat and diminfo in config
        rel_path = relpath(dirname(self._config_fn), self._mapdir)
        self.set_attribute_value('INPUT:CINPMAT', '"{:s}/inpmat-tmp.bin"'.format(rel_path))
        self.set_attribute_value('INPUT:LBMIROF', ".TRUE.")
        self.set_attribute_value('MAP:CDIMINFO', '"{:s}/diminfo_tmp.txt"'.format(rel_path))

    def write_config(self):
        """write adapted config to file. just before initializing
        only for models which do not allow for direct access to model config via bmi"""
        # rename old file if called 'input_flood.nam'
        glib.configcheck(self, self.logger)
        # write new file
        self._config_fn = join(dirname(self._config_fn), 'input_flood.nam')
        if isfile(self._config_fn):
            os.unlink(self._config_fn)
            self.logger.warn("{:s} file overwritten".format(self._config_fn))
        glib.configwrite(self._config, self._config_fn, encoding='utf-8')
        self.logger.info('Ini file written to {:s}'.format(self._config_fn))

# UTILS
def CMFtime_2_datetime(t):
    # internal CMF time is in minutes since 1850
    return datetime(1850, 1, 1) + timedelta(minutes = t)

class NamConfigParser(ConfigParser):
    def __init__(self, **kwargs):
        defaults = dict(comment_prefixes=('!', '/'),
                        inline_comment_prefixes=('!'),
                        delimiters=('='))
        defaults.update(**kwargs)
        super(NamConfigParser, self).__init__(**defaults)
        self.SECTCRE = re.compile(r"&N(?P<header>[^]]+)")

    def write(self, fp, space_around_delimiters=False):
        """Write an .ini-format representation of the configuration state.
        If `space_around_delimiters' is True (the default), delimiters
        between keys and values are surrounded by spaces.
        """
        super(NamConfigParser, self).write(fp, space_around_delimiters=space_around_delimiters)

    def _write_section(self, fp, section_name, section_items, delimiter):
        """Write a single section to the specified `fp'."""
        fp.write(u"&N{}\n".format(section_name))
        for key, value in section_items:
            if value.lower().strip('.') in ['true', 'false']:
                value = '.TRUE.' if value.lower().strip('.')=='true' else '.FALSE.'
            else:
                try:
                    float(value.replace('D', 'e'))
                except:
                    if not value.startswith('"'):
                        value = '"{}'.format(value)
                    if not value.endswith('"'):
                        value = '{}"'.format(value)
            value = self._interpolation.before_write(self, section_name, key, value)
            if value is not None or not self._allow_no_value:
                value = delimiter + str(value).replace('\n', '\n\t')
            else:
                value = ""
            fp.write("{}{}\n".format(key.upper(), value))
        fp.write("/\n")

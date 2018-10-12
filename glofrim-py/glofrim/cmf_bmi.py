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

from utils import setlogger, closelogger
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

    def __init__(self, engine, loglevel=logging.INFO, logger=None):
        self._bmi = _bmi(engine = engine)
        if logger:
            self.logger = logger.getChild(self._name)
        else:
            self.logger = setlogger(None, self._name, thelevel=loglevel)
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
        # reset model time to make sure it is consistent with the model
        self._dt = self.get_time_step()
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime

    def initialize(self, config_fn):
        if not hasattr(self, '_config'):
            self.initialize_config(config_fn)
        self.initialize_model()
            
    def update(self, dt=None):
        # dt in seconds. if not given model timestep is used
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        if (dt is not None) and (dt != self._dt.total_seconds()):
            _dt = timedelta(seconds=dt)
            if not glib.check_dts_divmod(_dt, self._dt):
                msg = "Invalid value for dt in comparison to model dt. Make sure a whole number of model timesteps ({}) fit in the given dt ({})"
                raise ValueError(msg.format(self._dt, _dt))
            self._bmi.update(dt)
            self._t += _dt
        else:
            self._bmi.update()
            self._t += self._dt
        self.logger.info('updated model to datetime {}'.format(self._t.strftime("%Y-%m-%d %H:%M")))

    def update_until(self, t, dt=None):
        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update(dt=dt)

    # not defined in CMF
    def spinup(self):
        """PCR specific spinup function"""
        raise NotImplementedError()

    def finalize(self):
        self.logger.info('finalize bmi. Close logger.')
        self._bmi.finalize()
        closelogger(self.logger)


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
        if self.initialized:
            dt = self._bmi.get_time_step()
        else:
            dt = int(self.get_attribute_value('CONF:DT'))
        self._dt = timedelta(**{self.get_time_units(): dt})
        return self._dt 
        
    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, fill_value=1e20, **kwargs):
        var = self._bmi.get_var(long_var_name, **kwargs)
        var = np.where(var == fill_value, np.nan, var)
        # return var with switched axis (fortran to python translation)
        return var.reshape(var.shape[::-1])

    def get_value_at_indices(self, long_var_name, inds, fill_value=1e20, **kwargs):
        return self.get_value(long_var_name, fill_value=fill_value, **kwargs).flat[inds]

    def set_value(self, long_var_name, src, **kwargs):
        grid_shape = self.get_grid().shape
        if (grid_shape[0] == src.shape[1]) and (grid_shape[1] == src.shape[0]):
            self.logger.info('src array has been flipped to fit CMF grid')
            # reverse grid (python to fortran translation)
            src = src.reshape(src.shape[::-1])
        src = np.where(np.isnan(src), 0, src).astype(self.get_var_type(long_var_name))
        self._bmi.set_var(long_var_name, src)

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        val = self.get_value(long_var_name)
        val.flat[inds] = src
        self._bmi.set_var(long_var_name, val)

    """
    Grid Information Functions
    """

    def get_grid(self, fn_params=r'params.txt', fn_locs=r'./hires/location.txt', **kwargs):
        if not hasattr(self, 'grid') or (self.grid is None):
            _root = dirname(self._config_fn)
            # read grid definition from params file
            fn_params = join(self._mapdir, fn_params)
            if not isfile(fn_params): raise IOError("Params.txt file not found {}".format(fn_params))
            self.logger.info('Getting Unit-Catchment Grid info based on {}'.format(basename(fn_params)))
            p = read_params(fn_params)
            # get hires index info file
            fn_locs = join(self._mapdir, fn_locs)
            if not isfile(fn_locs): raise IOError("locations.txt file for hires index not found {}".format(fn_locs))
            
            # get catchment unit grid instance 
            # TODO: set crs
            self.grid = UCGrid(p['transform'], p['height'], p['width'], fn_locs)
            
            # add drainage direction
            fn_nextxy = glib.getabspath(str(self.get_attribute_value('MAP:CNEXTXY').strip('"')), _root)
            if not isfile(fn_locs): raise IOError("nextxy file for drainage direction data not found {}".format(fn_locs))
            self.logger.info('Getting drainage direction from {}'.format(basename(fn_nextxy)))
            self.grid.set_dd(fn_nextxy, ddtype='nextxy', **kwargs)

            # set unit catchment mask based on drainage direction mask
            self.grid.set_mask(self.grid._dd.r.mask[0, :, :]==False)

        return self.grid

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, str):
            start_time = datetime.strptime(start_time, "%Y-%m-%d") 
        if not isinstance(start_time, datetime):
            raise ValueError("invalid date type")
        self._startTime = start_time
        self._t = start_time
        self.set_attribute_value('SIMTIME:ISYEAR', start_time.year)
        self.set_attribute_value('SIMTIME:ISMON', start_time.month)
        self.set_attribute_value('SIMTIME:ISDAY', start_time.day)
       
    def set_end_time(self, end_time):
        if isinstance(end_time, str):
            end_time = datetime.strptime(end_time, "%Y-%m-%d") 
        if not isinstance(end_time, datetime):
            raise ValueError("invalid date type")
        self._endTime = end_time
        self.set_attribute_value('SIMTIME:IEYEAR', end_time.year)
        self.set_attribute_value('SIMTIME:IEMON', end_time.month)
        self.set_attribute_value('SIMTIME:IEDAY', end_time.day)

    def set_out_dir(self, out_dir):
        self.set_attribute_value('OUTPUT:COUTDIR', abspath(out_dir))
        self._outdir = abspath(out_dir)

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

    def set_inpmat(self, bounds, res, olat='NtoS'):
        """Set the CMF inpmat file model based on the grid definition of upstream model"""
        if not isfile(join(self._mapdir, 'generate_inpmat')):
            raise ValueError('{} not found'.format(join(self._mapdir, 'generate_inpmat')))
        w, s, e, n = bounds
        # generate inpmat
        cmd = './generate_inpmat {} {} {} {} {} {:s}'
        cmd = cmd.format(abs(res), w, e, n, s, olat)
        # print(cmd)
        glib.subcall(cmd, cwd=self._mapdir)
        
    def set_inpmat_attrs(self):
        # set new inpmat and diminfo in config
        rel_path = relpath(self._mapdir, dirname(self._config_fn))
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

# CMF data utils
def read_params(fn, col_width=12):
    """parse params.txt file in CMF map folder to obtain low-res grid definition"""
    rename = {'grid number (north-south)': 'height',
              'grid number (east-west)': 'width',
              'west  edge [deg]': 'west',
              'north edge [deg]': 'north',
              'grid size  [deg]': 'res'}
    with open(fn, 'r') as txt:
        p = {line[col_width:].strip(): float(line[:col_width].strip())
                    for line in txt.readlines()}
    p = {rename[key]: p[key] for key in p if key in rename.keys()}
    p['transform'] = rasterio.transform.from_origin(p['west'], p['north'], p['res'], p['res'])
    return p
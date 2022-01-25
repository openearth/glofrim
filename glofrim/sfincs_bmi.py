import numpy as np
import sys
from configparser import ConfigParser
import logging
import os
from os.path import join, isfile, abspath, dirname, basename, relpath
from datetime import datetime, timedelta
from scipy.signal import convolve2d
import rasterio
import re

from bmi.wrapper import BMIWrapper as _bmi

from glofrim.utils import setlogger, closelogger
from glofrim.gbmi import GBmi
from glofrim.grids import RGrid
import glofrim.glofrim_lib as glib

class Sfincs(GBmi):
    """
    Glofrim implementation of the Sfincs BMI adaptor.
    """
    _name = 'Sfincs'
    _long_name = 'Sfincs'
    _version = ''
    _var_units = {'SGCQin': 'm3/s', 'dA': 'm2', 'H': 'm'}
    _input_var_names = ['SGCQin', 'dA']
    _output_var_names = ['SGCQin', 'H']
    _area_var_name = 'dA'
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
    def initialize_config(self, config_fn, config_defaults={}):
        if self.initialized:
            raise Warning("model already initialized, it's therefore not longer possible to initialize the config")
        # config settings
        defaults = {
            # 'tstart': '20000101 000000',
            # 'tref': '20000101 000000',
        }
        defaults.update(**config_defaults)
        self._config_fn = abspath(config_fn)
        self._config = glib.configread(self._config_fn, encoding='utf-8', cf=ParConfigParser(defaults=defaults))
        self._datefmt = "%Y%m%d %H%M%S"
        # model time
        self._dt = self.get_time_step()
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model paths
        _root = dirname(self._config_fn)
        self._mapdir = _root
        self._outdir = _root
        self.logger.info('Config initialized')

    def initialize_model(self):
        if not hasattr(self, '_config_fn'):
            raise Warning('Run initialize_config before initialize_model')
        # self.write_config() # write updated config to file as bmi does not allow direct access
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
            dt = timedelta(seconds=dt)
            # because of the adaptive timestep scheme do not check the dt value
            # if not glib.check_dts_divmod(dt, self._dt):
            #     msg = "Invalid value for dt in comparison to model dt. Make sure a whole number of model timesteps ({}) fit in the given dt ({})"
            #     raise ValueError(msg.format(self._dt, dt))
        else:
            dt = self._dt
        t_next = self.get_current_time() + dt
        i = 0
        while self._t < t_next:
            self._bmi.update()
            self._t = self.get_current_time()
            i += 1
        self.logger.info('updated model to datetime {} in {:d} iterations'.format(self._t.strftime("%Y-%m-%d %H:%M:%S"), i))

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
            ref_timestr = self.get_attribute_value('tref')
            refTime = datetime.strptime(ref_timestr, self._datefmt)
            startTime_float = self._bmi.get_start_time()
            startTime = refTime + timedelta(**{self.get_time_units(): startTime_float})
        else:
            start_timestr = self.get_attribute_value('tstart')
            startTime = datetime.strptime(start_timestr, self._datefmt)
        # startTime = refdate + timedelta(**{self.get_time_units(): TStart})
        self._startTime = startTime
        return self._startTime
    
    def get_current_time(self):
        if self.initialized:
            curtime = timedelta(**{self.get_time_units(): self._bmi.get_current_time()})
            return self._startTime + curtime
        else:
            return self.get_start_time()

    def get_end_time(self):
        if self.initialized:
            ref_timestr = self.get_attribute_value('tref')
            refTime = datetime.strptime(ref_timestr, self._datefmt)
            endTime_float = self._bmi.get_end_time()
            endTime = refTime + timedelta(**{self.get_time_units(): endTime_float})

        else:
            stop_timestr = self.get_attribute_value('tstop')
            endTime = datetime.strptime(stop_timestr, self._datefmt)

        self._endTime = endTime
        return self._endTime

    def get_time_step(self):
        if self.initialized:
            dt = self._bmi.get_time_step()
        else:
            dt = float(1.e-6)  # TODO: check if this is
        self._dt = timedelta(**{self.get_time_units(): dt})
        return self._dt 

    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    def get_mask(self, long_var_name):
        if long_var_name in ['Qx', 'QxSGold', 'Qy', 'QySGold']:
            if long_var_name in ['Qx', 'QxSGold']:
                # expand mask in x-direction
                mask = np.vstack([self.grid.mask, np.zeros((1, self.grid.mask.shape[1]))])
                conv_filt = np.array([[1, 1]])  # horizontal convolution
            elif long_var_name in ['Qy', 'QySGold']:
                mask = np.hstack([self.grid.mask, np.zeros((self.grid.mask.shape[0], 1))])
                conv_filt = np.array([[1], [1]])  # vertical convolution
            mask = np.array(convolve2d(mask, conv_filt), dtype='bool')
        else:
            mask = self.grid.mask
        return mask

    def get_value(self, long_var_name, **kwargs):
        var = np.rot90(np.asarray(self._bmi.get_var(long_var_name)).copy())
        mask = self.get_mask(long_var_name)
        var[~mask] = np.nan
        # ensure boundary vals are also set to NaN
        var[var==-99999] = np.nan
        return var

    def get_value_at_indices(self, long_var_name, inds, **kwargs):
        return self.get_value(long_var_name).flat[inds]

    def set_value(self, long_var_name, src, fill_value=0., **kwargs):
        # set nans that lie within to_mod model domain to zeros to prevent model crashes
        mask = self.get_mask(long_var_name)
        src[mask & np.isnan(src)] = 0.
        # set remaining nans to missing value
        src = np.where(np.isnan(src), fill_value, src).astype(self.get_var_type(long_var_name))
        # LFP does not have a set_var function, but used the get_var function with an extra argument
        self._bmi.set_var(long_var_name, np.rot90(src, k=3))  # rotate to match the order of Sfincs internal grids

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        val = self.get_value(long_var_name)
        val.flat[inds] = src
        self.set_value(long_var_name, val)

    """
    Grid Information Functions
    """
    def get_grid(self):
        if not hasattr(self, 'grid') or (self.grid is None):
            # get the 2d dem values directly from bmi
            dem_grid = np.rot90(self._bmi.get_var("zb"))
            transform = rasterio.transform.Affine(
                float(self.get_attribute_value('dx')),
                0.,
                float(self.get_attribute_value('x0')),
                0.,
                float(self.get_attribute_value('dy')),
                float(self.get_attribute_value('y0')),
            )
            self.grid = RGrid(
                transform,
                int(self.get_attribute_value('nmax')) + 2,  # the +2 is for the boundary conditions
                int(self.get_attribute_value('mmax')) + 2,  # the +2 is for the boundary conditions
                crs=int(self.get_attribute_value('epsg')),
                mask=np.isfinite(dem_grid)
            )
            # # riv width file used for "1D coords"
            # _width_fn = glib.getabspath(str(self.get_attribute_value('SGCwidth')), self._mapdir)
            # if not isfile(_width_fn): raise IOError('SGCwidth file not found')
            # with rasterio.open(_width_fn, 'r') as ds:
            #     row, col = np.where(ds.read(1)>0)
            #     x, y = self.grid.xy(row=row, col=col)
            #     inds = self.grid.ravel_multi_index(row, col)
            # self.grid.set_1d(nodes=np.vstack((x, y)).transpose(), links=None, inds=inds)  # python2.7 nodes=np.array(zip(x, y))
        return self.grid


    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, datetime):
            refdate = start_time.strftime(self._datefmt)
        elif isinstance(start_time, str):
            try:
                refdate = start_time # str
                start_time = datetime.strptime(start_time, self._datefmt) # check format
            except ValueError:
                raise ValueError('wrong date format, use "yyyy-mm-dd"')
        else:
            raise ValueError('wrong start_date datatype')
        self._startTime = start_time
        self._t = start_time
        self.set_attribute_value('refdate', refdate)

    def set_end_time(self, end_time):
        if isinstance(end_time, str):
            try:
                end_time = datetime.strptime(end_time, self._datefmt)
            except ValueError:
                raise ValueError('wrong end_date format, use "yyyy-mm-dd"')
        if not isinstance(end_time, datetime):
            raise ValueError('wrong end_date datatype')
        refdate = self.get_start_time()
        assert end_time >  refdate
        TStop = (end_time - refdate).seconds + (end_time - refdate).days * 86400
        TStop = '{:.0f}'.format(TStop)
        self._endTime = end_time
        self.set_attribute_value('sim_time', TStop)

    def set_out_dir(self, out_dir):
        self.set_attribute_value('dirroot', relpath(out_dir, dirname(self._config_fn)))
        self._outdir = abspath(out_dir)

    def get_attribute_names(self):
        glib.configcheck(self, self.logger)
        return glib.configattr(self._config)
    
    def get_attribute_value(self, attribute_name):
        glib.configcheck(self, self.logger)
        # always use "general" as config header; as file has no config header this is hard-coded
        if ':' not in attribute_name:
            attribute_name = 'general:{}'.format(attribute_name)
        else:
            attribute_name = 'general:{}'.format(attribute_name.split(':')[1])
        self.logger.debug("get_attribute_value: {}".format(attribute_name))
        return glib.configget(self._config, attribute_name)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        glib.configcheck(self, self.logger)
        # always use "general" as config header; as file has no config header this is hard-coded
        if ':' not in attribute_name:
            attribute_name = 'general:{}'.format(attribute_name)
        else:
            attribute_name = 'general:{}'.format(attribute_name.split(':')[1])
        self.logger.debug("set_attribute_value: {} -> {}".format(attribute_name, attribute_value))
        return glib.configset(self._config, attribute_name, str(attribute_value))

    def write_config(self):
        """write adapted config to file. just before initializing
        only for models which do not allow for direct access to model config via bmi"""
        self._config_fn = glib.write_config(self, self._config, self._config_fn, self.logger)

# UTILS
class ParConfigParser(ConfigParser):
    def __init__(self, **kwargs):
        self.optionxform = lambda option:option # keep format with capital/lower letters
        defaults = dict(comment_prefixes=('!', '/', '#'),
                        inline_comment_prefixes=('!'), allow_no_value=True,
                        delimiters=('='))
        defaults.update(**kwargs)
        super(ParConfigParser, self).__init__(**defaults)

    def read_file(self, f, **kwargs):
        def par2ini(f, header_name):
            """change par to ini before parse as ini
            note that this removes comments"""
            yield '[{}]\n'.format(header_name)
            for line in f:
                # yield '='.join(line.split()[:2])
                yield line
        super(ParConfigParser, self).read_file(par2ini(f, 'general'), **kwargs)
        
    def _write_section(self, fp, section_name, section_items, delimiter):
        """Write a single section to the specified `fp'."""
        for key, value in section_items:
            value = self._interpolation.before_write(self, section_name, key, value)
            value = ' ' + str(value).replace('\n', '\n\t')
            fp.write("{}{}\n".format(key, value))
        fp.write("\n")

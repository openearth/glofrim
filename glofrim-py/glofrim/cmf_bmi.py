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
from gbmi import GBmi, GBmiGridType, GBmiModelType
import gbmi_lib as glib 

class CMF(GBmi):
    """
    Glofrim implementation of the PCR BMI adaptor.
    """
    _name = 'CaMa-Flood'
    _version = '3.6.2'
    _var_units = {'roffin': 'mm', # NOTE: unit set by user
                  'outflw': 'm3/s', }
    _input_var_names = ['roffin']
    _output_var_names = ['outflw']
    _area_var_name = ''

    def __init__(self, engine):
        self._bmi = _bmi(engine = engine)
        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False

    """
    Model Control Functions
    """
    def initialize_config(self, config_fn):
        if self.initialized:
            raise Warning("model already initialized, it's therefore not longer possible to initialize the config")
        # config settings
        self._config_fn = abspath(config_fn)
        self._config = glib.configread(self._config_fn, encoding='utf-8', cf=NamConfigParser())
        # model time
        self._dt = timedelta(seconds=int(self.get_attribute_value('CONF:DT')))
        self._timeunit = 'seconds'
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model files
        root = dirname(self._config_fn)
        self._nextxy_fn = glib.getabspath(str(self.get_attribute_value('MAP:CNEXTXY').strip('"')), root)
        self._nextxy_fn = self._nextxy_fn.replace('.bin', '.tif')
        self._indir = dirname(self._nextxy_fn)
        self._outdir = glib.getabspath(str(self.get_attribute_value('OUTPUT:COUTDIR').strip('"')), root)
        if not isfile(self._nextxy_fn): raise IOError('nextxy file {} not found'.format(self._nextxy_fn))
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
    Model Information Functions
    """
    def get_model_type(self):
        return GBmiModelType.ROUTING 

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
        return self._dt 
        
    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, **kwargs):
        return np.asarray(self._bmi.get_var(long_var_name))

    def get_value_at_indices(self, long_var_name, inds, **kwargs):
        # always use 1d inds
        if not ((isinstance(inds, np.ndarray)) and (inds.ndim == 1)):
            raise ValueError('indices should be 1d arrays')
        return self.get_value(long_var_name).flat[inds]

    def set_value(self, long_var_name, src, **kwargs):
        self._bmi.set_var(long_var_name, src)

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        # always use 1d inds
        if not ((isinstance(inds, np.ndarray)) and (inds.ndim == 1)):
            raise ValueError('indices should be 1d arrays')
        val = self.get_value(long_var_name)
        val.flat[inds] = src
        self._bmi.set_var(long_var_name, val)

    def get_drainage_direction(self, **kwargs):
        from nb.nb_io import read_dd_rasterio
        self.logger.info('Getting drainage direction from {}'.format(basename(self._nextxy_fn)))
        self.dd = read_dd_rasterio(self._nextxy_fn, ddtype='nextxy', **kwargs)

    """
    Grid Information Functions
    """
    def get_grid_type(self):
        return GBmiGridType.UNITCATCHMENT

    def get_grid_transform(self):
        if not hasattr(self, 'grid_transform'):
            self.logger.info('Getting grid transform based on {}'.format(basename(self._nextxy_fn)))
            with rasterio.open(self._nextxy_fn, 'r') as ds:
                self.grid_transform = ds.transform
        return self.grid_transform

    def get_grid_shape(self):
        if not hasattr(self, 'grid_shape'):
            self.logger.info('Getting grid shape based on {}'.format(basename(self._nextxy_fn)))
            with rasterio.open(self._nextxy_fn, 'r') as ds:
                self.grid_shape = ds.shape
        return self.grid_shape

    def get_grid_bounds(self):
        if not hasattr(self, 'grid_bounds'):
            self.logger.info('Getting grid bounds based on {}'.format(basename(self._nextxy_fn)))
            with rasterio.open(self._nextxy_fn, 'r') as ds:
                self.grid_bounds = ds.bounds
        return self.grid_bounds

    def get_grid_res(self):
        if not hasattr(self, 'grid_res'):
            self.logger.info('Getting grid bounds based on {}'.format(basename(self._nextxy_fn)))
            with rasterio.open(self._nextxy_fn, 'r') as ds:
                self.grid_res = ds.res
        return self.grid_res

    def grid_index(self, x, y, catmx_fn=r'./hires/reg.catmxy.tif', **kwargs):
        """Get CMF indices (1d) of xy coordinates using the catmxy.tif file"""
        fn_catmxy = glib.getabspath(catmx_fn, self._indir)
        if not isfile(fn_catmxy):
            raise IOError("{} file not found".format(fn_catmxy))
        # read catmxy temporary into memory
        self.logger.info('Getting grid index based on {}'.format(basename(fn_catmxy)))
        with rasterio.open(fn_catmxy, 'r', driver='GTiff') as ds:
            ncount = ds.count
            if ncount != 2:
                raise ValueError("{} file should have two layers".format(fn_catmxy))
            # python zero based index for high res CMF grid
            rows, cols = ds.index(x, y)
            rows, cols = np.atleast_1d(rows).astype(int), np.atleast_1d(cols).astype(int)
            # make sure indices are inside hr grid
            nrows_hr, ncols_hr = ds.shape
            inside_hr = np.logical_and.reduce((rows>=0,rows<nrows_hr,cols>=0,cols<ncols_hr))
            c, r = np.ones_like(cols)*-1, np.ones_like(rows)*-1
            # read low-res CMF fortran one-based index
            c[inside_hr], r[inside_hr] = ds.read()[:, rows[inside_hr, None], cols[inside_hr, None]].squeeze()
            # go from fortran one-based to python zero-based indices
            r[inside_hr], c[inside_hr] = r[inside_hr]-1, c[inside_hr]-1
        # check if inside domain
        nrows, ncols = self.get_grid_shape()
        inside = np.logical_and.reduce((r>=0, r<nrows, c>=0, c<ncols))
        # calculate 1d index 
        # NOTE invalid indices have value -1
        inds = np.ones_like(inside, dtype=int)*-1 
        if np.any(inside):
            inds[inside] = np.ravel_multi_index((r[inside], c[inside]), (nrows, ncols))
        return inds

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, str):
            try:
                start_time.strptime("%Y-%m-%d") 
            except ValueError:
                raise ValueError('wrong date format, use "yyyy-mm-dd"')
        if not isinstance(start_time, datetime):
            raise ValueError("invalid date type")
        self.set_attribute_value('SIMTIME:ISYEAR', start_time.year)
        self.set_attribute_value('SIMTIME:ISMON', start_time.month)
        self.set_attribute_value('SIMTIME:ISDAY', start_time.day)
       
    def set_end_time(self, end_time):
        if isinstance(end_time, str):
            try:
                end_time.strptime("%Y-%m-%d") 
            except ValueError:
                raise ValueError('wrong date format, use "yyyy-mm-dd"')
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

    def write_config(self):
        """write adapted config to file. just before initializing
        only for models which do not allow for direct access to model config via bmi"""
        # rename old file if called 'input_flood.nam'
        glib.configcheck(self, self.logger)
        if basename(self._config_fn) == 'input_flood.nam':
            os.rename(self._config_fn, self._config_fn.replace('.nam', '.temp'))
        # write new file
        self._config_fn = join(dirname(self._config_fn), 'input_flood.nam')
        if isfile(self._config_fn):
            os.unlink(self._config_fn)
            self.logger.warn("{:s} file overwritten".format(self._config_fn))
        glib.configwrite(self._config, self._config_fn, encoding='utf-8')
        self.logger.info('Ini file written to {:s}'.format(self._config_fn))

    def set_inpmat(self, bounds, res, olat='NtoS', DROFUNIT=1):
        """Set the CMF inpmat file model based on the grid definition of upstream model"""
        self.logger.info("Setting inpmat file")
        ddir = self._indir
        if not isfile(join(ddir, 'generate_inpmat')):
            raise ValueError('{} not found'.format(join(ddir, 'generate_inpmat')))
        if not abs(res[0]) == abs(res[1]):
            raise ValueError('lat and lon resolution should be the same in regular grid')
        westin  = bounds.left
        eastin  = bounds.right
        northin = bounds.top
        southin = bounds.bottom
        # generate inpmat
        cmd = './generate_inpmat {} {} {} {} {} {:s}'
        cmd = cmd.format(abs(res[0]), westin, eastin, northin, southin, olat)
        self.logger.info(cmd)
        glib.subcall(cmd, cwd=ddir)
        # set new inpmat and diminfo in config
        rel_path = relpath(dirname(self._config_fn), ddir)
        self.set_attribute_value('INPUT:CINPMAT', '"{:s}/inpmat-tmp.bin"'.format(rel_path))
        self.set_attribute_value('INPUT:LBMIROF', ".TRUE.")
        self.set_attribute_value('MAP:CDIMINFO', '"{:s}/diminfo_tmp.txt"'.format(rel_path))
        self.set_attribute_value('CONF:DROFUNIT', '{:d}'.format(int(DROFUNIT)))

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

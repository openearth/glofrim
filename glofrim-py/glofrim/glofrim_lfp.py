# -*- coding: utf-8 -*-

import glob
import shutil
import os
from os import mkdir
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs
import numpy as np
import rtree
import StringIO
import rasterio
from configparser import ConfigParser
import re
from bmi.wrapper import BMIWrapper
from main import BMI_model_wrapper

class LFP_model(BMI_model_wrapper):
    def __init__(self, engine, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date,
                 missing_value=np.nan, **kwargs):
        """initialize the LISFLOOD-FP (LFP) model BMI class and model configuration file"""
        # TODO: extend this list to cover all variables
        si_unit_conversions = {'rain': 1e-3, ## [mm/s] -> [m/s]
                               } ##
        ## initialize BMIWrapper and model
        lfp_bmi = BMIWrapper(engine = engine)

        # set config parser        
        configparser = ParConfigParser()

        # for offline use the forcing data dir can be set. not yet inplemented
        forcing_data_dir = ''
        options = dict(dt=86400, tscale=1., # sec / dt
                        missing_value=missing_value)
        # initialize BMIWrapper for model
        super(LFP_model, self).__init__(lfp_bmi, config_fn, 'LFP', 'sec',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, configparser=configparser,
                                        si_unit_conversions=si_unit_conversions,
                                        **kwargs)

        # set some basic model properties
        duration = (end_date - start_date).total_seconds()
        
        #TODO: needs to be changed to fit with structure of par-files
#        globalOptions = {'dummysection':
#							{"sim_time": "{:f}".format(duration),},
#						}
        globalOptions = {}
        
        self.set_config(globalOptions)

    def initialize(self):
        # move model input files to out dir
        self.set_model_input_files()
        # write updated config and intialize
        super(LFP_model, self).initialize()
        # read mesh after initialization
        self.get_model_coords()

    def set_model_input_files(self):
        src = self.model_data_dir
        dst = self.out_dir
        for fn in glob.glob(src + '/*'):
            if isfile(fn):
                shutil.copy(fn, dst)
            elif isdir(fn):
                if not isdir(join(dst, basename(fn))):
                    mkdir(join(dst, basename(fn)))

    ## model grid functionality
    def get_model_coords(self):
        """Get LFP model coordinates for 1D and 2D mesh via BMI. The LFP model
        should be initialized first in order to access the variables."""

        self.logger.info('Getting LFP model coordinates.')

        i_ind, j_ind = np.where(np.logical_and(self.get_var('SGCwidth') > 0., self.get_var('DEM') != -9999))
        # print i_ind.shape, j_ind.shape
        
        fn_map = join(self.model_data_dir,
                      self.model_config['header1']['DEMfile'])

        if not isfile(fn_map):
            raise IOError('landmask file not found')
        self._fn_landmask = fn_map

        with rasterio.open(fn_map, 'r') as ds:
            self.grid_index = ds.index
            self.model_grid_res = ds.res
            self.model_grid_bounds = ds.bounds
            self.model_grid_shape = ds.shape
            self.model_grid_transform = ds.transform
            list_x_coords, list_y_coords = ds.xy(i_ind, j_ind)

        self.model_1d_coords = zip(list_x_coords, list_y_coords)
        self.model_1d_indices = np.arange(i_ind.size)
        self.model_1d_rc = (i_ind, j_ind)

        pass

    def get_area_1d(self):
        row, col = self.model_1d_rc
        area_1D = self.get_var('dA')[row, col]
        return area_1D


class ParConfigParser(ConfigParser):
    def __init__(self, **kwargs):
        self.optionxform = lambda option:option # keep format with capital/lower letters
        defaults = dict(comment_prefixes=('!', '/', '#'),
                        inline_comment_prefixes=('!'),
                        delimiters=('='))
        defaults.update(**kwargs)
        super(ParConfigParser, self).__init__(**defaults)

    def read_file(self, f, **kwargs):
        def par2ini(f, header_name):
            """change par to ini before parse as ini
            note that this removes comments"""
            yield '[{}]\n'.format(header_name)
            for line in f:
                yield '='.join(line.split()[:2])
        super(ParConfigParser, self).read_file(par2ini(f, 'header1'), **kwargs)
        
    def _write_section(self, fp, section_name, section_items, delimiter):
        """Write a single section to the specified `fp'."""
        for key, value in section_items:
            value = self._interpolation.before_write(self, section_name, key, value)
            value = ' ' + str(value).replace('\n', '\n\t')
            fp.write("{}{}\n".format(key, value))
        fp.write("\n")
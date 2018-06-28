# -*- coding: utf-8 -*-

import logging
import glob
import shutil
from os import mkdir
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs
import numpy as np
import rtree

from bmi.wrapper import BMIWrapper

from main import BMI_model_wrapper
from utils import ConfigParser

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)


class DFM_model(BMI_model_wrapper):
    def __init__(self, engine, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date,
                 missing_value=np.nan, **kwargs):
        """initialize the Delft3D-FM (DFM) model BMI class and model configuration file"""
        # TODO: extend this list to cover all variables
        si_unit_conversions = {'rain': 1e-3, ## [mm/s] -> [m/s]
                               } ##
        ## initialize BMIWrapper and model
        dfm_bmi = BMIWrapper(engine = engine)
        # set config parser
        configparser = ConfigParser(inline_comment_prefixes=('#'))
        # for offline use the forcing data dir can be set. not yet inplemented
        forcing_data_dir = ''
        options = dict(dt=86400, tscale=1., # sec / dt
                        missing_value=missing_value)
        # initialize BMIWrapper for model
        super(DFM_model, self).__init__(dfm_bmi, config_fn, 'DFM', 'sec',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, configparser=configparser,
                                        si_unit_conversions=si_unit_conversions,
                                        **kwargs)
        # set some basic model properties
        globalOptions = {'time':
                            {'RefDate': start_date.strftime("%Y%m%d"),
                             'TStart': 0,
                             'TStop': int((end_date - start_date).total_seconds())
                            },
                         'output':
                            {'OutputDir': ""} # use default output dir settings
                        }
        self.set_config(globalOptions)


    def initialize(self):
        # move model input files to out dir
        self.set_model_input_files()
        # write updated config and intialize
        super(DFM_model, self).initialize()
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
        """Get DFM model coordinates for 1D and 2D mesh via BMI. The DFM model
        should be initialized first in order to access the variables."""

        logger.info('Getting DFM model coordinates.')
        # define separator between 2D and 1D parts of arrays == lenght of 2d cell points
        self._1d2d_idx = len(self.get_var('flowelemnode'))
        x_coords = self.get_var('xz') # x-coords of each cell centre point
        y_coords = self.get_var('yz') # y-coords of each cell centre point
        xy_coords = zip(x_coords, y_coords)
        self.model_2d_coords = xy_coords[:self._1d2d_idx]
        self.model_2d_indices = range(self._1d2d_idx)
        self.model_1d_coords = xy_coords[self._1d2d_idx:]
        n1d = len(self.model_1d_coords)
        self.model_1d_indices = np.arange(n1d, dtype=np.int32) + self._1d2d_idx
        pass

    def get_area_1d(self):
        #NOTE: to stay consistent with the indices we should return the full array
        return self.get_var('ba')

    def get_model_1d_index(self):
        """Creat a spatial index for the 1d coordinates. A model_1d_index
        attribute funtion is created to find the nearest 1d coordinate tuple"""
        logger.info('Constructing spatial index for the 1D coordinates of the DFM model.')
        # 1d coords
        n1d = len(self.model_1d_coords)
        self.model_1d_indices = np.arange(n1d, dtype=np.int32) + self._1d2d_idx
        # build spatial rtree index of points2
        logger.info('Constructing spatial index for 1D vertices of DFM model')
        self.model_1d_rtree = rtree.index.Index()
        for i, xy in enumerate(self.model_1d_coords):
            self.model_1d_rtree.insert(i+n1d, xy) # return index including 2d

    def model_1d_index(self, xy, n=1):
        if not hasattr(self, 'model_1d_rtree'):
            self.get_model_1d_index()
        xy = [xy] if isinstance(xy, tuple) else xy
        return [list(self.model_1d_rtree.nearest(xy0, 1))[0] for xy0 in xy]

    def get_model_2d_index(self):
        """Creat a spatial index for the 2d mesh center coordinates.
        A model_2d_index attribute funtion is created to find the nearest
        2d cell center"""
        # build spatial rtree index of 2d coords
        logger.info('Constructing spatial index for the 2D mesh of the DFM model')
        self.model_2d_rtree = rtree.index.Index()
        for i, xy in enumerate(self.model_2d_coords):
            self.model_2d_rtree.insert(i, xy)

    def model_2d_index(self, xy, n=1):
        if not hasattr(self, 'model_2d_rtree'):
            self.get_model_2d_index()
        xy = [xy] if isinstance(xy, tuple) else xy
        index = [list(self.model_2d_rtree.nearest(xy0, 1))[0] for xy0 in xy]
        #TODO: assess validity based on e.g. max distance
        valid = np.ones(len(index), dtype=bool)
        return index, valid

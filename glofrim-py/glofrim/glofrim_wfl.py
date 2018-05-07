# -*- coding: utf-8 -*-

import logging
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs
import numpy as np
import rasterio

# from pcrglobwb_bmi_v203 import pcrglobwb_bmi
from wflow import wflow_bmi

from main import BMI_model_wrapper
from utils import ConfigParser

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)

class WFL_model(WFL_model_wrapper):
    def __init__(self, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date,
                 missing_value=-999, landmask_mv=255, forcing_data_dir=None,
                 **kwargs):
        """initialize the WFLOW (WFL) model BMI class and model configuration file"""
        # BMIWrapper for PCR model model
        pcr_bmi = pcrglobwb_bmi.pcrglobwbBMI()
        # set config parser
        configparser = ConfigParser(inline_comment_prefixes=('#'))
        # model and forcing data both in model_data_dir
        if forcing_data_dir is None:
            forcing_data_dir = model_data_dir
        options = dict(dt=1, tscale=86400., # seconds per dt
                        missing_value=missing_value, landmask_mv=landmask_mv)
        # initialize BMIWrapper for model
        super(PCR_model, self).__init__(pcr_bmi, config_fn, 'PCR', 'day',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, configparser=configparser, **kwargs)
        # set some basic model properties
        globalOptions = {'globalOptions':
                            {'inputDir': self.forcing_data_dir,
                             'outputDir': self.out_dir,
                             'startTime': start_date.strftime("%Y-%m-%d"),
                             'endTime': end_date.strftime("%Y-%m-%d")
                        }}
        self.set_config(globalOptions)

        # read grid and ldd properties
        self.get_model_grid()
        self.get_drainage_direction()

    def get_model_grid(self):
        """Get WFL model bounding box, resolution and index based on the
        landmask map file.

        Created attributes
        -------
        model_grid_res : tuple
            model grid x, y resolution
        model_grid_bounds : list
            model grid xmin, ymin, xmax, ymax bounds
        model_grid_shape : tuple
            model number of rows and cols
        """
        logger.info('Getting WFLOW model grid parameters.')
        fn_map = join(self.model_config['globalOptions']['inputDir'],
                      self.model_config['globalOptions']['landmask'])
        if not isfile(fn_map):
            raise IOError('landmask file not found')
        self._fn_landmask = fn_map
        with rasterio.open(fn_map, 'r') as ds:
            self.grid_index = ds.index
            self.model_grid_res = ds.res
            self.model_grid_bounds = ds.bounds
            self.model_grid_shape = ds.shape
            self.model_grid_transform = ds.transform
        msg = 'Model bounds {:s}; width {}, height {}'
        logger.debug(msg.format(self.model_grid_bounds, *self.model_grid_shape))
        pass

    def get_drainage_direction(self):
        from nb.nb_io import read_dd_pcraster
        # read file with pcr readmap
        nodata = self.options.get('landmask_mv', 255)
        logger.info('Getting PCR LDD map.')
        fn_ldd = self.model_config['routingOptions']['lddMap']
        if not isabs(fn_ldd):
            ddir = self.model_config['globalOptions']['inputDir']
            fn_ldd = abspath(join(ddir, fn_ldd))
        if not isfile(fn_ldd):
            raise IOError('ldd map file {} not found.'.format(fn_ldd))
        self.dd = read_dd_pcraster(fn_ldd, self.model_grid_transform, nodata=nodata)
        pass

    def model_2d_index(self, xy, **kwargs):
        """Get PCR (row, col) indices of xy coordinates.

        Arguments
        ----------
        xy : list of tuples
          list of (x, y) coordinate tuples

        Returns
        -------
        indices : list of tuples
          list of (row, col) index tuples
        """
        import pcraster as pcr
        nodata = self.options.get('landmask_mv', 255)
        logger.info('Getting PCR model indices of xy coordinates.')
        r, c = self.grid_index(*zip(*xy), **kwargs)
        r = np.array(r).astype(int)
        c = np.array(c).astype(int)
        # get valid cells (inside landmask)
        lm_map = pcr.readmap(str(self._fn_landmask))
        lm_data = pcr.pcr2numpy(lm_map, nodata)
        valid = lm_data[r, c] == 1
        return zip(r, c), valid

    def get_coupled_flux(self):
        """Get summed runoff and upstream discharge at coupled cells"""
        if not hasattr(self, 'coupled_mask'):
            msg = 'The PCR model must be coupled before the total coupled flux can be calculated'
            raise AssertionError(msg)
        # get PCR runoff and discharge at masked cells
        runoff = np.where(self.coupled_mask == 1, np.nan_to_num(self.get_var('runoff')) * self.get_var('cellArea'), 0) # [m3/day]
        q_out = np.where(self.coupled_mask == 2,  np.nan_to_num(self.get_var('discharge')) * 86400., 0) # [m3/day]
        # sum runoff + discharge routed one cell downstream
        tot_flux = runoff + self.dd.dd_flux(q_out)
        return tot_flux * self.options['dt']

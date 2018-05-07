# -*- coding: utf-8 -*-

#TODO: guess the input dir as in PCR has to be set up to work with path of the WFLOW ini-file instead

import logging
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs
import numpy as np
import rasterio

# from pcrglobwb_bmi_v203 import pcrglobwb_bmi
import wflow.bmi
from wflow.wflow_bmi import wflowbmi_light, wflowbmi_csdms

from main import BMI_model_wrapper
from utils import ConfigParser

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)

class WFL_model(BMI_model_wrapper):
    def __init__(self, config_fn,
                 out_dir,
                 start_date, end_date,
                 missing_value=-999, landmask_mv=255, forcing_data_dir=None,
                 **kwargs):
        """initialize the WFLOW (WFL) model BMI class and model configuration file"""

        # BMIWrapper for WFLOW model model
        model_data_dir = dirname(config_fn)

        # set config parser
        configparser = ConfigParser(inline_comment_prefixes=('#'))
        
        # model and forcing data both in model_data_dir
        if forcing_data_dir is None:
            forcing_data_dir = model_data_dir
        options = dict(dt=1, tscale=86400., # seconds per dt
                        missing_value=missing_value, landmask_mv=landmask_mv)
        # initialize BMIWrapper for model
        super(WFL_model, self).__init__(wflowbmi_csdms(), config_fn, 'WFL', 'day',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, configparser=configparser, **kwargs)
        # set some basic model properties
        # TODO: set model start and end time
        # globalOptions = {'run': {'startTime': start_date.strftime("%Y-%m-%d"),
        #                          'endTime': end_date.strftime("%Y-%m-%d")
        #                 }}
        # self.set_config(globalOptions)

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
        fn_map = getattr(self.model_config['model'], 'wflow_subcatch', 'staticmaps/wflow_subcatch.map')
        if not isabs(fn_map):
            ddir = self.model_data_dir
            fn_map = abspath(join(ddir, fn_map))
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
        logger.info('Getting WFLOW LDD map.')
        fn_ldd = getattr(self.model_config['model'], 'wflow_ldd', 'staticmaps/wflow_ldd.map')
        if not isabs(fn_ldd):
            ddir = self.model_data_dir
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
        logger.info('Getting WFLOW model indices of xy coordinates.')
        r, c = self.grid_index(*zip(*xy), **kwargs)
        r = np.array(r).astype(int)
        c = np.array(c).astype(int)
        # get valid cells (inside landmask)
        with rasterio.open(self._fn_landmask, 'r') as src:
            lm_data = np.ma.masked_equal(src.read(1), src.nodata)
            valid = lm_data.mask[r, c] == False
        return zip(r, c), valid

    def get_coupled_flux(self):
        """Get summed runoff and upstream discharge at coupled cells"""
        if not hasattr(self, 'coupled_mask'):
            msg = 'The WFL model must be coupled before the total coupled flux can be calculated'
            raise AssertionError(msg)
        # get PCR runoff and discharge at masked cells
        runoff = np.where(self.coupled_mask == 1, np.nan_to_num(self.get_var('runoff')) * self.get_var('cellArea'), 0) # [m3/day]
        q_out = np.where(self.coupled_mask == 2,  np.nan_to_num(self.get_var('discharge')) * 86400., 0) # [m3/day]
        # sum runoff + discharge routed one cell downstream
        tot_flux = runoff + self.dd.dd_flux(q_out)
        return tot_flux * self.options['dt']

    ## exchange states
    # NOTE: overwrite parent bmi get_var en set_var because naming is different
    def get_var(self, name, parse_missings=True, mv=None):
        var = self.bmi.get_attribute_value(name)
        if parse_missings: # if given nodata is parsed to np.nan
            mv = self._mv if mv is None else mv
            var = np.where(var == mv, np.nan, var)
        if name in self._si_unit_conversions: # convert model var to SI units
            var = var * float(self._si_unit_conversions[name])
        return var

    def set_var(self, name, var, parse_missings=True, mv=None):
        if name in self._si_unit_conversions: # convert var from SI to model units
            var = var / float(self._si_unit_conversions[name])
        if parse_missings: # set nans back to model mv data values
            mv = self._mv if mv is None else mv
            var = np.where(np.isnan(var), mv, var)
        self.bmi.set_attribute_value(name, var)

    # NOTE: overwrite parent bmi update because it takes not dt argument
    def update(self, **kwargs):
        """Advance model state by one time step.

        Perform all tasks that take place within one pass through the model's
        time loop. This typically includes incrementing all of the model's
        state variables.
        """
        self.bmi.update()
        current_time = self.get_current_time()
        time_step = self.get_time_step()
        logger.info(
            "%s -> start_time: %s, current_time %s, timestep %s",
            self.name,
            self.start_time,
            current_time,
            time_step
        )
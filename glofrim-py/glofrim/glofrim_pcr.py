# -*- coding: utf-8 -*-

import logging
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs
import numpy as np
import rasterio

from pcrglobwb_bmi_v203 import pcrglobwb_bmi

from main import BMI_model_wrapper
from utils import ConfigParser

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)

class PCR_model(BMI_model_wrapper):
    def __init__(self, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date,
                 missing_value=-999, landmask_mv=255, forcing_data_dir=None,
                 **kwargs):
        """initialize the PCR-GLOBWB (PCR) model BMI class and model configuration file"""
        # BMIWrapper for PCR model model
        pcr_bmi = pcrglobwb_bmi.pcrglobwbBMI()
        # set config parser
        self._configparser = ConfigParser(inline_comment_prefixes=('#'))
        # model and forcing data both in model_data_dir
        if forcing_data_dir is None:
            forcing_data_dir = model_data_dir
        options = dict(dt=1, tscale=86400., # seconds per dt
                        missing_value=missing_value, landmask_mv=landmask_mv)
        # initialize BMIWrapper for model
        super(PCR_model, self).__init__(pcr_bmi, config_fn, 'PCR', 'day',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, **kwargs)
        # set some basic model properties
        globalOptions = {'globalOptions':
                            {'inputDir': self.forcing_data_dir,
                             'outputDir': self.out_dir,
                             'startTime': start_date.strftime("%Y-%m-%d"),
                             'endTime': end_date.strftime("%Y-%m-%d")
                        }}
        self.set_config(globalOptions)

    def get_model_grid(self):
        """Get PCR model bounding box, resolution and index based on the
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
        from nb.nb_io import read_dd_pcraster
        logger.info('Getting PCR model grid parameters.')
        fn_map = join(self.model_config['globalOptions']['inputDir'],
                      self.model_config['globalOptions']['landmask'])
        if not isfile(fn_map):
            raise IOError('landmask file not found')
        with rasterio.open(fn_map, 'r') as ds:
            self._model_index = ds.index
            self.model_grid_res = ds.res
            self.model_grid_bounds = ds.bounds
            self.model_grid_shape = ds.shape
            self.model_grid_transform = ds.transform
        self._fn_landmask = fn_map
        msg = 'Model bounds {:s}; width {}, height {}'
        logger.debug(msg.format(self.model_grid_bounds, *self.model_grid_shape))

        # read file with pcr readmap
        logger.info('Getting PCR LDD map.')
        fn_ldd = self.model_config['routingOptions']['lddMap']
        if not isabs(fn_ldd):
            ddir = self.model_config['globalOptions']['inputDir']
            fn_ldd = abspath(join(ddir, fn_ldd))
        if not isfile(fn_ldd):
            raise IOError('ldd map file {} not found.'.format(fn_ldd))
        nodata = self.options.get('landmask_mv', 255)
        self.dd = read_dd_pcraster(fn_ldd, self.model_grid_transform, nodata=nodata)

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
        logger.info('Getting PCR model indices of xy coordinates.')
        if not hasattr(self, '_model_index'):
            self.get_model_grid()
        r, c = self._model_index(*zip(*xy), **kwargs)
        r = np.array(r).astype(int)
        c = np.array(c).astype(int)
        # get valid cells (inside landmask)
        nodata = self.options.get('landmask_mv', 255)
        lm_map = pcr.readmap(str(self._fn_landmask))
        lm_data = pcr.pcr2numpy(lm_map, nodata)
        valid = lm_data[r, c] == 1
        return zip(r, c), valid


    def deactivate_routing(self, coupled_indices=None):
        """Deactive LDD at indices by replacing the local ldd value with 5, the
        ldd pit value. The ldd is modified offline. This only has effect before
        the model is initialized.

        Arguments
        ----------
        coupled_indices : list of tuples, str
          list with (x, y) tuples or 'all' to deactivate total grid
        """
        # this function requires pcr functionality
        import pcraster as pcr
        if coupled_indices is None:
            if not hasattr(self, 'coupled_dict'):
                msg = 'The PCR model must be coupled before deactivating ' + \
                      'routing in the coupled cells'
                raise AssertionError(msg)
            coupled_indices = self.coupled_dict.keys()
        if self.initialized:
            msg = "Deactivating the LDD is only possible before the model is initialized"
            raise AssertionError(msg)

        logger.info('Editing PCR ldd grid to deactivate routing in coupled cells.')
        # get ldd filename from config
        fn_ldd = self.model_config['routingOptions']['lddMap']
        if not isabs(fn_ldd):
            ddir = self.model_config['globalOptions']['inputDir']
            fn_ldd = abspath(join(ddir, fn_ldd))
        # read file with pcr readmap
        if not isfile(fn_ldd):
            raise IOError('ldd map file {} not found.'.format(fn_ldd))
        lddmap = pcr.readmap(str(fn_ldd))
        # change nextxy coupled indices to pits == 5
        nodata = self.options.get('landmask_mv', 255)
        ldd = pcr.pcr2numpy(lddmap, nodata)
        if coupled_indices == 'all':
            ldd = np.where(ldd != nodata, 5, ldd)
        else:
            rows, cols = zip(*coupled_indices)
            rows, cols = np.atleast_1d(rows), np.atleast_1d(cols)
            ldd[rows, cols] = np.where(ldd[rows, cols] != nodata, 5, ldd[rows, cols])
        # write to tmp file
        self._fn_ldd_tmp = join(self.out_dir, basename(fn_ldd))
        # write map with pcr.report
        lddmap = pcr.numpy2pcr(pcr.Ldd, ldd, nodata)
        pcr.report(lddmap, str(self._fn_ldd_tmp))
        # set tmp file in config
        ldd_options = {'routingOptions': {'lddMap': self._fn_ldd_tmp}}
        self.set_config(ldd_options)

    def get_coupled_flux(self):
        """Get summed runoff and upstream discharge at coupled cells"""
        if not hasattr(self, 'coupled_mask'):
            msg = 'The PCR model must be coupled before the total coupled flux can be calculated'
            raise AssertionError(msg)
        # get PCR runoff and discharge at masked cells
        runoff = np.where(self.coupled_mask == 1,  self.get_var('runoff') * self.get_var('cellArea'), 0) # [m3/day] 
        q_out = np.where(self.coupled_mask == 2,  self.get_var('discharge') * 86400., 0) # [m3/day] 
        # sum runoff + discharge routed one cell downstream
        tot_flux = runoff + self.dd.dd_flux(q_out)
        return tot_flux * self.options['dt']


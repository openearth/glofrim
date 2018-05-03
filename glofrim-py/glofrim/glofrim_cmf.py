# -*- coding: utf-8 -*-

import logging
import glob
import shutil
from os import mkdir, unlink
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs, splitext, relpath
from datetime import datetime, timedelta
import numpy as np
import rasterio

from bmi.wrapper import BMIWrapper

from main import BMI_model_wrapper
from utils import NamConfigParser, subcall

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)


class CMF_model(BMI_model_wrapper):
    def __init__(self, engine, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date,
                 missing_value=1e20, **kwargs):
        """initialize the CaMa-Flood (CMF) model BMI class and model configuration file"""
        ## initialize BMIWrapper and model
        cmf_bmi = BMIWrapper(engine = engine)
        # set config parser
        self._configparser = NamConfigParser()
        # for offline use the forcing data dir can be set. not yet inplemented
        forcing_data_dir = ''
        options = dict(dt=86400, tscale=1, # sec / dt
                        missing_value=missing_value)
        # initialize BMIWrapper for model
        super(CMF_model, self).__init__(cmf_bmi, config_fn, 'CMF', 'sec',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, **kwargs)
        # setup output dir
        if not isdir(join(self.out_dir, 'out')):
            mkdir(join(self.out_dir, 'out'))
        # set some basic model properties
        globalOptions = {'NSIMTIME':
                            {"ISYEAR": "{:d}".format(start_date.year),
                            "ISMON": "{:d}".format(start_date.month),
                            "ISDAY": "{:d}".format(start_date.day),
                            "IEYEAR": "{:d}".format(end_date.year),
                            "IEMON": "{:d}".format(end_date.month),
                            "IEDAY": "{:d}".format(end_date.day)
                            },
                        'NOUTPUT':
                            {'COUTDIR': '"{}/out/"'.format(self.out_dir)},
                        }
        self.set_config(globalOptions)

    def initialize(self):
        # move model input files to out dir
        self.set_model_input_files()
        # write updated config and intialize
        super(CMF_model, self).initialize()

    def set_model_input_files(self):
        # move files
        move_dict = {'NMAP': 'map', 'NINPUT': 'input'}
        mdir = dirname(self.config_fn)
        for sec in move_dict:
            folder = move_dict[sec]
            dst_path = join(self.out_dir, folder)
            if not isdir(dst_path):
                mkdir(dst_path)
            for opt in self.model_config[sec]:
                fpath = self.model_config[sec][opt].strip('"')
                fn = basename(fpath)
                if not isabs(fpath): # path relative to config_fn
                    rel_path = dirname(fpath)
                    fpath = join(mdir, rel_path, fn)
                if isfile(fpath):
                    # copy all files with same name, ignore extensions
                    for src_fn in glob.glob('{}.*'.format(splitext(fpath)[0])):
                        shutil.copy(src_fn, dst_path)
                        # remove *tmp.* files from data_dir
                        if 'tmp.' in  basename(src_fn):
                            logger.info('removing tmp file {:s} from model data dir'.format(basename(src_fn)))
                            unlink(src_fn)
                    # update config_fn
                    self.update_config(sec, opt, '"./{}/{}"'.format(folder, fn))

    def set_inpmat_file(self, bounds, res, olat='NtoS'):
        """Set the CMF inpmat file model based on the grid definition of upstream
        model"""
        ddir = self.model_data_dir
        if not isfile(join(ddir, 'generate_inpmat')):
            raise ValueError('{} not found'.format(join(ddir, 'generate_inpmat')))
        if not abs(res[0]) == abs(res[1]):
            raise ValueError('lat and lon resolution should be the same in regular grid')
        westin  = bounds.left
        eastin  = bounds.right
        northin = bounds.top
        southin = bounds.bottom
        # generate inpmat
        msg2 = './generate_inpmat {} {} {} {} {} {:s}'.format(
                        abs(res[0]), westin, eastin, northin, southin, olat)
        logger.info(msg2)
        subcall(msg2, cwd=ddir)
        # set new inpmat and diminfo in config
        rel_path = relpath(dirname(self.config_fn), ddir)
        inpmat_options = {'NINPUT': {'CINPMAT': '"{:s}/inpmat-tmp.bin"'.format(rel_path),
                                     'LBMIROF': '.TRUE.'
                                    },
                        'NMAP': {'CDIMINFO': '"{:s}/diminfo_tmp.txt"'.format(rel_path)
                                },
                        'NCONF': {'DROFUNIT': '1'} #  SI units [m]
                        }
        self.set_config(inpmat_options)

    ## model grid functionality
    def get_model_grid(self):
        """Get CMF model bounding box, resolution and shape based on landmask
        geotiff file 'lsmask.tif'.

        To convert cama bin maps to geotiff use the 'cama_maps_io.py' script.

        The function creates the following attributes
        -------
        model_grid_res : tuple
            model grid x, y resolution
        model_grid_bounds : list
            model grid xmin, ymin, xmax, ymax bounds
        model_grid_shape : tuple
            model number of rows and cols
        """
        from nb.dd_ops import NextXY

        logger.info('Getting CMF model grid parameters.')
        fn_lsmask = join(self.model_data_dir, 'lsmask.tif')
        if not isfile(fn_lsmask):
            raise IOError("lsmask.tif file not found at {}".format(fn_lsmask))
        # set model low-res grid parameters
        with rasterio.open(fn_lsmask, 'r') as ds:
            self.model_grid_res = ds.res
            self.model_grid_bounds = ds.bounds
            self.model_grid_shape = ds.shape
            self.model_grid_transform = ds.transform
        self._fn_landmask = fn_lsmask
        msg = 'Model bounds {:s}; width {}, height {}'
        logger.debug(msg.format(self.model_grid_bounds, *self.model_grid_shape))

        # read drainage direction data
        logger.info('Getting CMF model drainage direction')
        fn_nextxy = join(self.model_data_dir, 'nextxy.bin')
        if not isfile(fn_nextxy):
            raise IOError("nextxy.bin file not found at {}".format(fn_nextxy))
        # read drainage direction data and initialize nextxy object.
        nextxy_data = np.fromfile(fn_nextxy, dtype=np.int32).reshape(2, *self.model_grid_shape)
        self.dd = NextXY(nextxy_data, transform=self.model_grid_transform, nodata=self._mv)

    def model_2d_index(self, xy, **kwargs):
        """Get CMF (row, col) indices at low resolution based on xy coordinates.
        The low-res indices are looked up catmxy geotiff file 'reg.catmxy.tif'.

        To convert cama binary maps to geotiff use the 'cama_maps_io.py' script.

        Note that CMF indices smaller than zero should be ignored! This occurs
        for unitcatchments that are not represented in CMF or fall out of
        the CMF domain. for both (row, col) = (-1, -1). CMF is volume conservative,
        so the runoff of ignored unitcatchments is conserved by spreading it
        of other unitcatchments in the same cell.

        Arguments
        ---------
        xy : list of tuples
          list of (x, y) coordinate tuples

        Returns
        -------
        indices : list of tuples
          list of (row, col) index tuples
        """
        logger.info('Getting CMF model indices for xy coordinates.')
        fn_catmxy = join(self.model_data_dir, 'hires', 'reg.catmxy.tif')
        if not isfile(fn_catmxy):
            raise IOError("{} file not found".format(fn_catmxy))
        # read catmxy temporary into memory
        with rasterio.open(fn_catmxy, 'r', driver='GTiff') as ds:
            ncount = ds.count
            if ncount != 2:
                raise ValueError("{} file should have two layers".format(fn_catmxy))
            # python zero based index for high res CMF grid
            rows, cols = ds.index(*zip(*xy))
            rows, cols = np.atleast_1d(rows).astype(int), np.atleast_1d(cols).astype(int)
            # make sure indices are inside grid
            nrows_hr, ncols_hr = ds.shape
            inside = np.logical_and.reduce((rows>=0,rows<nrows_hr,cols>=0,cols<ncols_hr))
            cmf_cols, cmf_rows = np.zeros_like(cols), np.zeros_like(rows)
            # read low-res CMF fortran one-based index
            cmf_cols[inside], cmf_rows[inside] = ds.read()[:, rows[inside, None], cols[inside, None]].squeeze()
            # go from fortran one-based to python zero-based indices
            cmf_rows, cmf_cols = cmf_rows-1, cmf_cols-1
        # check valid indices -> only cells within domain // row, col -1 values should be ignored
        nrows, ncols = self.model_grid_shape
        valid = np.logical_and.reduce((cmf_rows>=0,cmf_rows<nrows,cmf_cols>=0,cmf_cols<ncols))
        return zip(cmf_rows, cmf_cols), valid
        
    def get_var(self, name, parse_missings=True, *args, **kwargs):
        var = super(CMF_model, self).get_var(name, parse_missings=parse_missings)
        # return var with switched axis (fortran to python translation)
        return var.reshape(var.shape[::-1])

    def get_current_time(self):
        t = super(CMF_model, self).get_current_time()
        return self._CMFtime_2_datetime(t)

    def get_start_time(self):
        t = super(CMF_model, self).get_start_time()
        return self._CMFtime_2_datetime(t)

    def _CMFtime_2_datetime(self, t):
        # internal CMF time is in minutes since 1850
        return datetime(1850, 1, 1) + timedelta(minutes = t)

    def get_coupled_flux(self):
        """Get summed runoff and upstream discharge at coupled cells"""
        if not hasattr(self, 'coupled_mask'):
            msg = 'The CMF model must be coupled before the total coupled flux can be calculated'
            raise AssertionError(msg)
        # get CMF runoff and discharge at masked cells
        runoff = np.where(self.coupled_mask == 1, np.nan_to_num(self.get_var('runoff')), 0) # [m3/s]
        q_out = np.where(self.coupled_mask == 2,  np.nan_to_num(self.get_var('outflw')), 0) # [m3/s]
        # sum runoff + discharge routed one cell downstream
        tot_flux = runoff + self.dd.dd_flux(q_out)
        return tot_flux * self.options['dt']
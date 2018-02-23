# -*- coding: utf-8 -*-
# TODO: write main docstring
# TODO: add decorator to child functions to get docstring from parent function

import rasterio
import numpy as np
from pcrglobwb_bmi_v203 import pcrglobwb_bmi
from bmi.wrapper import BMIWrapper
import logging
import warnings
from datetime import datetime, timedelta
from collections import OrderedDict
import os, shutil
import re
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs
import glob
import rtree

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)

# local libraries
from utils import subcall, config_to_dict, dict_to_config, ConfigParser, NamConfigParser

# wrapper around BMI
class BMI_model_wrapper(object):
    def __init__(self, bmi, config_fn, name, t_unit,
                 model_data_dir, forcing_data_dir, out_dir,
                 options, si_unit_conversions={}, **kwargs):
        """initialize the model BMI class and model configuration file"""
        # TODO: extend docstring
        self.bmi = bmi
        self.name = name # string used
        self.t_unit = t_unit
        # set model paths
        self.config_fn = abspath(config_fn)
        self.model_data_dir = abspath(model_data_dir)
        self.forcing_data_dir = abspath(forcing_data_dir)
        self.out_dir = abspath(out_dir)
        if not isdir(self.out_dir):
            os.mkdir(self.out_dir)
        # first step of two step initialization.
        self.initialize_config()
        # second step is not yet performed
        self.initialized = False
        # set some class specific options and internal shortcuts
        self.options = options
        self._dt = options['dt']
        self._mv = options['missing_value']
        # convert all units to SI (lenght -> meters, time -> seconds)
        self._si_unit_conversions = si_unit_conversions


    # model configuration
    def initialize_config(self, config_fn=None, **kwargs):
        """Read text based model configuration file to internal model_config
        dictionary.

        First step of two-phase initialization. In this step only the configuration
        is read in. This allows a user to then change settings and parameters
        before fully initializing the model

        Arguments
        ----------
        config_fn : str, optional
          The path to the model configuration file. If given it overwrites the
          internal config_fn attribute.
        """
        if config_fn is not None:
            self.config_fn = abspath(config_fn)
        self.model_config = config_to_dict(self.config_fn,
                                           cf=self._configparser,
                                           **kwargs)

    def write_config(self, **kwargs):
        """The internal model_config dictionary is written to the out_dir. This
        step should be excecuted just before the model initialization."""
        self.config_fn = join(self.out_dir, basename(self.config_fn))
        dict_to_config(self.model_config, self.config_fn,
                       cf=self._configparser, **kwargs)
        logger.info('Ini file for {:s} written to {:s}'.format(self.name, self.config_fn))

    def set_config(self, model_config):
        """Change multiple model config file settings with dictionary.
        This will only have affect before the model is Initialized.

        Arguments
        ----------
        model_config : dictionary
          A nested dictionary with section, setting and value to be replaced.

        Example input
        -------------
        {SECTION1:
            {setting1: value1,
             setting2: value2},
        SECTION2:
            {setting3: value3}
            }
        """

        if not isinstance(model_config, dict):
            raise ValueError("input not of dictionary type")
        else:
            if not isinstance(model_config, dict):
                msg = "input requires nested dictionary type, see function documentation"
                raise ValueError(msg)
        # update internal config dictionary
        for sec in model_config:
            for opt in model_config[sec]:
                value = model_config[sec][opt]
                self.update_config(sec, opt, value)

    def update_config(self, sec, opt, value):
        """Change model config file settings. This will only have affect before
        the model is initialized.

        Arguments
        ----------
        sec : str
          Configuration file section header name
        opt : str
          Configuration file option name
        value : str
          Configuration file option value
        """
        self.model_config[sec].update(**{opt: value})

    # set class options
    # TODO: not sure this adds much and makes it less transparent. Should perhaps
    # be replaced at some point
    def set_options(self, **kwargs):
        self.options.update({kw: arg for kw in kwargs})
        if 'dt' in kwargs:
            self._dt = self.options['dt']
        if 'missing_value' in kwargs:
            self._mv = self.options['missing_value']

    # model initialization
    def initialize(self):
        """Perform startup tasks for the model.

        This is second step of two-phase initialization and includes writing the
        internal model configuration dictionary to file and reading the model
        grid or mesh coordinates.
        """
        # write possibly updated config file
        self.write_config()
        # initialize model with updated config file
        self.bmi.initialize(self.config_fn)
        # set start time attribute
        self.start_time = self.get_start_time()
        self.initialized = True
        logger.info('{:s} initialized'.format(self.name))

    def finalize(self):
        """Shutdown the library and clean up the model.
        Note that the out_dir is not cleaned up.
        """
        self.bmi.finalize()

    def spinup(self):
        """Spin-up the model in order to create realistic initiale state.
        Note that this is not implemented for every model."""
        if not hasattr(self.bmi, 'spinup'):
            raise NotImplementedError("Spin-up functionality not implemented")
        self.bmi.spinup(*args, **kwargs)

    # exchange states
    def get_var(self, name, parse_missings=True, mv=None):
        """Get data as nd array from the given variable <name>. All variables are
        returned in SI units based on the internal _si_unit_conversions dict.

        Arguments
        ----------
        name : str
          An input or output variable name. use the get_var_name_all function to
          get a list with all exposed model variable names.
        parse_missings : bool (default True)
          If True, missing values are parsed to np.nan based on internal missing
          value attribute
        mv : int, float, optional
          Use temporal value to parse missing values

        Returns
        -------
        nd array
          numpy array with data
        """
        var = self.bmi.get_var(name)
        if parse_missings: # if given nodata is parsed to np.nan
            mv = self._mv if mv is None else mv
            var = np.where(var == mv, np.nan, var)
        if name in self._si_unit_conversions: # convert model var to SI units
            var = var * float(self._si_unit_conversions[name])
        return var

    def set_var(self, name, var, parse_missings=True, mv=None):
        """Write nd array data to given variable <name> model state. All variables
        should be set in SI units and are converted to model units based on the
        internal _si_unit_conversions dict.

        Arguments
        ----------
        name : str
          An input or output variable name. use the get_var_name_all function to
          get a list with all exposed model variable names.
        var : nd array
            numpy array with data
        parse_missings : bool (default True)
          If True, np.nan values are parsed to the model missing value based on
          the internal missing value attribute
        mv : int, float, optional
          Use temporal value to parse np.nan values to missing values
        """

        if name in self._si_unit_conversions: # convert var from SI to model units
            var = var / float(self._si_unit_conversions[name])
        if parse_missings: # set nans back to model mv data values
            mv = self._mv if mv is None else mv
            var = np.where(np.isnan(var), mv, var)
        self.bmi.set_var(name, var)

    def set_var_index(self, name, index, var, parse_missings=True, mv=None):
        """Write nd array data at specific indices to given variable <name> model
        state. All variables should be set in SI units and are converted to model
        units based on the internal _si_unit_conversions dict.

        Arguments
        ----------
        name : str
          An input or output variable name. use the get_var_name_all function to
          get a list with all exposed model variable names.
        index : list, tuple
          For 1D input: a list with indices
          For 2D input: list of (row, col) tuples OR a tuple of row and col lists
        var : nd array
          Numpy nd array with varialbe data
        parse_missings : bool (default True)
          If True, np.nan values are parsed to the model missing value based on
          the internal missing value attribute
        mv : int, float, optional
          Use temporal value to parse np.nan values to missing values
        """

        if isinstance(index[0], tuple):
            index = zip(*index) # from list of x, y tuples to tuple of x, y lists
        if name in self._si_unit_conversions: # convert var from SI to model units
            var = var / float(self._si_unit_conversions[name])
        if parse_missings: # set nans back to model mv data values
            mv = self._mv if mv is None else mv
            var = np.where(np.isnan(var), mv, var)
        # bug in bmi wrapper line 639
        # if hasattr(self.bmi, 'set_var_index'):
        #     self.bmi.set_var_index(name, index, var)
        # else: # PCR has no set_var_index function
        var0 = self.get_var(name, parse_missings=False)
        var0[index] = var
        self.bmi.set_var(name, var0)

    # run timestep dt in model
    def update(self, dt=None):
        """Advance model state by one time step.

        Perform all tasks that take place within one pass through the model's
        time loop. This typically includes incrementing all of the model's
        state variables.
        """
        if dt is None: # by default take internally set dt
            dt = self._dt
        self.bmi.update(dt=dt)
        current_time = self.get_current_time()
        time_step = self.get_time_step()
        logger.info(
            "%s -> start_time: %s, current_time %s, timestep %s",
            self.name,
            self.start_time,
            current_time,
            time_step
        )

    # time attributes
    def get_start_time(self):
        "get model start (initialization) time"
        return self.bmi.get_start_time()

    def get_current_time(self):
        "get model current time"
        return self.bmi.get_current_time()

    def get_time_step(self):
        "get model current time step"
        return self.bmi.get_time_step()

    # var info
    def get_var_count(self):
        return self.bmi.get_var_count()

    def get_var_name(self, i):
        return self.bmi.get_var_name(i)

    def get_var_name_all(self):
        return [self.get_var_name(i) for i in xrange(self.get_var_count())]

    ## coupling functions
    def couple_grid_to_1d(self, other):
        """Couple external 1d coordinates to internal model 2d grid. The model
        routing is deactivated at coupled cells and area fractions are calculated
        to distribute water volume accross 1d nodes.

        The output model grid indices are model specific:
        - PCR : regular grid (row, col)
        - CMF : irregular unit catchment grid (row, col)
        - DFM : flexible mesh (flat index)

        Arguments
        ---------
        other : BMI_model_wrapper object
            downstream BMI_model_wrapper object

        Created Attributes
        ------------------
        coupled_idx : list
            All model indices of the coupled cells (upstream model) to 1d nodes (downstream model)
        coupled_dict : dict
            Per coupled model index, the indices of the other model
        coupled_area_frac : nd array
            Fraction of water volume from upstream model grid cell that should be
            added to downstream 1d node
        """
        if (other.name != 'DFM') or (self.name not in ['PCR', 'CMF']):
            msg = 'Grid to 1D coupling has only been implemented for PCR to' + \
                  ' CMF (upstream) to DFM (downstream) coupling'
            raise NotImplementedError(msg)

        if other.name == 'DFM':
            if not other.initialized:
                msg = 'The DFM model should be initialized first to obtain ' + \
                      ' the 1D model coordinates via BMI'
                raise AssertionError(msg)
            if not hasattr(other, 'model_1d_coords'):
                other.get_model_coords()
            xy_other, indices_other = other.model_1d_coords, other.model_1d_indices
            area_other = other.get_var('ba')

        logger.info('Coupling {} grid to {} 1D nodes.'.format(self.name, other.name))
        # get cell indices at 1D coordinates
        cellidx = self.model_2d_index(xy_other)
        # set indices to easily exchange variables
        other.coupled_idx = indices_other
        self.coupled_idx = zip(*cellidx) # tuple of (row, col) lists
        # create coupled 1-to-1 downstream-to-upstream indices dictionary
        other.coupled_dict = {i: idx for i, idx in zip(indices_other, cellidx)}
        # invert dictionary for 1-to-n upstream-to-downstream coupling
        self.coupled_dict = dictinvert(other.coupled_dict)

        logger.info('Getting fraction of coupled 1d nodes based on area.')
        # get area fractions of coupled cells
        area_frac = {} # [m2/m2]
        for _, idx in list(self.coupled_dict.items()):
            area_other_sum_idx = np.sum(area_other[idx])
            afs = {i: area_other[i]/area_other_sum_idx for i in idx}
            area_frac.update(afs)
        self.coupled_area_frac = np.array([area_frac[i] for i in indices_other])

        # deactivate routing at coupled cells
        self.deactivate_routing()

    def couple_grid_to_grid(self, other):
        """Couple external grid to internal model 2d grid. The model
        ldd grid is deactivated at coupled cells.

        The exact grid coupling method is model specific:
        - CMF: via the inpmat
        Note: only implemented for the CMF model

        Arguments
        ---------
        other : BMI_model_wrapper object
            downstream BMI_model_wrapper object
        """
        if (self.name != 'PCR') or (other.name != 'CMF'):
            msg = 'Grid to grid coupling has only been implemented PCR to CMF (other) coupling'
            raise NotImplementedError(msg)

        # check model grid extends
        if not hasattr(self, 'model_grid_bounds'):
            self.get_model_grid()
        if not hasattr(other, 'model_grid_bounds'):
            other.get_model_grid()
        bounds, res = self.model_grid_bounds, self.model_grid_res

        logger.info('Coupling {} grid to {} grid.'.format(self.name, other.name))
        # set impmat of CMF model based on the model grid bounds and resolution
        if other.name == 'CMF':
            other.set_inpmat_file(bounds, res)

        # deactivate routing in upstream model
        self.deactivate_routing('all')


class PCR_model(BMI_model_wrapper):
    def __init__(self, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date, dt=1,
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
        options = dict(dt=dt, tscale=86400, # seconds per dt
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
        self._fn_landmask = fn_map
        msg = 'Model bounds {:s}; width {}, height {}'
        logger.debug(msg.format(self.model_grid_bounds, *self.model_grid_shape))

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

        logger.info('Getting PCR model indices of xy coordinates.')
        if not hasattr(self, '_model_index'):
            self.get_model_grid()
        r, c = self._model_index(*zip(*xy), **kwargs)
        r = np.array(r).astype(int)
        c = np.array(c).astype(int)
        return zip(r, c)

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

class CMF_model(BMI_model_wrapper):
    def __init__(self, engine, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date, dt=86400,
                 missing_value=1e20, **kwargs):
        """initialize the CaMa-Flood (CMF) model BMI class and model configuration file"""
        ## initialize BMIWrapper and model
        cmf_bmi = BMIWrapper(engine = engine)
        # set config parser
        self._configparser = NamConfigParser()
        # for offline use the forcing data dir can be set. not yet inplemented
        forcing_data_dir = ''
        options = dict(dt=dt, tscale=1, # sec / dt
                        missing_value=missing_value)
        # initialize BMIWrapper for model
        super(CMF_model, self).__init__(cmf_bmi, config_fn, 'CMF', 'sec',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, **kwargs)
        # setup output dir
        if not isdir(join(self.out_dir, 'out')):
            os.mkdir(join(self.out_dir, 'out'))
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
                            {'COUTDIR': '"./out"'},
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
                os.mkdir(dst_path)
            for opt in self.model_config[sec]:
                fpath = self.model_config[sec][opt].strip('"')
                fn = basename(fpath)
                if not os.path.isabs(fpath): # path relative to config_fn
                    rel_path = dirname(fpath)
                    fpath = join(mdir, rel_path, fn)
                if isfile(fpath):
                    # copy all files with same name, ignore extensions
                    for src_fn in glob.glob('{}.*'.format(os.path.splitext(fpath)[0])):
                        shutil.copy(src_fn, dst_path)
                        # remove *tmp.* files from data_dir
                        if 'tmp.' in  basename(src_fn):
                            logger.info('removing tmp file {:s} from model data dir'.format(basename(src_fn)))
                            os.unlink(src_fn)
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
        rel_path = os.path.relpath(dirname(self.config_fn), ddir)
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

        logger.info('Getting CMF model grid parameters.')
        fn_lsmask = join(self.model_data_dir, 'lsmask.tif')
        if not isfile(fn_lsmask):
            raise IOError("lsmask.tif file not found at {}".format(fn_lsmask))
        # set model low-res grid parameters
        with rasterio.open(fn_lsmask, 'r') as ds:
            self.model_grid_res = ds.res
            self.model_grid_bounds = ds.bounds
            self.model_grid_shape = ds.shape
        self._fn_landmask = fn_lsmask
        msg = 'Model bounds {:s}; width {}, height {}'
        logger.debug(msg.format(self.model_grid_bounds, *self.model_grid_shape))

    def model_2d_index(self, xy, **kwargs):
        """Get CMF (row, col) indices at low resolution based on xy coordinates.
        The low-res indices are looked up catmxy geotiff file 'reg.catmxy.tif'.

        To convert cama binary maps to geotiff use the 'cama_maps_io.py' script.

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
                raise ValueError("{} file should have two layers".format(fn))
            nrows, ncols = ds.shape
            rows, cols = ds.index(*zip(*xy))
            # go from fortran one-based to python zero-based indices
            rows, cols = np.atleast_1d(rows).astype(int)-1, np.atleast_1d(cols).astype(int)-1
            invalid = np.logical_and.reduce((rows<0,rows>=nrows,cols<0,cols>=ncols))
            if np.any(invalid):
                raise IndexError('XY coordinates outside of CMF domain')
            cmf_idx = ds.read()[:, rows[:, None], cols[:, None]].squeeze()
            cmf_cols, cmf_rows = cmf_idx
        return zip(cmf_rows, cmf_cols)

    def deactivate_routing(self, coupled_indices=None):
        """Deactive nextxy at coupled indices by replacing the nextxy index
        values with -9, the nextxy pit value.

        Arguments
        ----------
        coupled_indices : list of tuples, str
          list with (row, col) tuples of low-res CMF model
        """
        if coupled_indices is None:
            if not hasattr(self, 'coupled_dict'):
                msg = 'The CMF model must be coupled before deactivating ' + \
                      'routing in the coupled cells'
                raise AssertionError(msg)
            coupled_indices = self.coupled_dict.keys()
        if self.initialized:
            msg = "Deactivate routing is only possible before the model is initialized"
            raise AssertionError(msg)
        # model_grid_shape attribute required to read data correctly
        if not hasattr(self, 'model_grid_shape'):
            self.get_model_grid()

        logger.info('Editing CMF nextxy grid to deactivate routing in coupled cells.')
        # read input data
        fn_nextxy = join(self.model_data_dir, 'nextxy.bin')
        if not isfile(fn_nextxy):
            raise IOError("lsmask.tif file not found at {}".format(fn_nextxy))
        nextxy = np.fromfile(fn_nextxy, dtype=np.int32).reshape(2, *self.model_grid_shape)
        rows, cols = zip(*coupled_indices)
        # change nextxy coupled indices to pits == -9
        rows, cols = np.atleast_1d(rows), np.atleast_1d(cols)
        nextxy[:, rows[:, None], cols[:, None]] = -9
        # get file name for tmp file based on orignal file name from config
        fn_nextxy_org = self.model_config['NMAP']['CNEXTXY'].strip('"')
        if not isabs(fn_nextxy_org):
            fn_nextxy_org = abspath(join(dirname(self.config_fn), fn_nextxy_org))
        fn_nextxy_new = fn_nextxy_org.replace('.bin', '_tmp.bin')
        # write data to tmp file
        nextxy.astype(np.int32).tofile(fn_nextxy_new)
        # set new nextxy in config
        rel_path = os.path.relpath(dirname(self.config_fn), dirname(fn_nextxy_new))
        nextxy_options = {'NMAP': {'CNEXTXY': join(rel_path, basename(fn_nextxy_new))}}
        self.set_config(nextxy_options)

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

class DFM_model(BMI_model_wrapper):
    def __init__(self, engine, config_fn,
                 model_data_dir, out_dir,
                 start_date, end_date, dt=86400,
                 missing_value=np.nan, **kwargs):
        """initialize the Delft3D-FM (DFM) model BMI class and model configuration file"""
        # TODO: extend this list to cover all variables
        si_unit_conversions = {'rain': 1e-3, ## [mm/s] -> [m/s]
                               } ##
        ## initialize BMIWrapper and model
        dfm_bmi = BMIWrapper(engine = engine)
        # set config parser
        self._configparser = ConfigParser(inline_comment_prefixes=('#'))
        # for offline use the forcing data dir can be set. not yet inplemented
        forcing_data_dir = ''
        options = dict(dt=dt, tscale=1., # sec / dt
                        missing_value=missing_value)
        # initialize BMIWrapper for model
        super(DFM_model, self).__init__(dfm_bmi, config_fn, 'DFM', 'sec',
                                        model_data_dir, forcing_data_dir, out_dir,
                                        options, si_unit_conversions=si_unit_conversions,
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

    def get_model_1d_index(self):
        """Creat a spatial index for the 1d coordinates. A model_1d_index
        attribute funtion is created to find the nearest 1d coordinate tuple"""

        if not hasattr(self, 'model_1d_coords'):
            self.get_model_coords()
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
        if not hasattr(self, 'model_2d_coords'):
            self.get_model_coords()
        # build spatial rtree index of 2d coords
        logger.info('Constructing spatial index for the 2D mesh of the DFM model')
        self.model_2d_rtree = rtree.index.Index()
        for i, xy in enumerate(self.model_2d_coords):
            self.model_2d_rtree.insert(i, xy)

    def model_2d_index(self, xy, n=1):
        if not hasattr(self, 'model_2d_rtree'):
            self.get_model_2d_index()
        xy = [xy] if isinstance(xy, tuple) else xy
        return [list(self.model_2d_rtree.nearest(xy0, 1))[0] for xy0 in xy]


# utils
def dictinvert(d):
    inv = {}
    for k, v in d.iteritems():
        keys = inv.setdefault(v, [])
        keys.append(k)
    return inv

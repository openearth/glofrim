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
from os.path import isdir, join, basename, dirname, abspath, isfile
import glob
import rtree

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)

# local libraries
from utils import subcall, config_to_dict, dict_to_config, ConfigParser, NamConfigParser

# wrapper around BMI
class _model(object):
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

        Parameters
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

        Parameters
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

        Parameters
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

        Parameters
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

        Parameters
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

        Parameters
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

    def couple_1d_2d(self, xy, indices=None):
        """Couple external 1d coordinates to internal model 2d grid. A dictionary
        with for each 1d xy coordinates, the index of the 2d grid (1 to 1) and
        its inversed dictionary (1 to n).

        The 2d grid is dependent of the model and therefore the internal 2d indices
        - PCR : regular grid (row, col)
        - CMF : irregular unit catchment grid (row, col)
        - DFM : flexible mesh (flat index)

        Parameters
        ----------
        xy : list of tuples
          list of (x, y) coordinate tuples
        indices : list or nd array, optional
          if provided these indices are used to create the output dictionary
        """
        if not hasattr(self, 'model_2d_index'):
            logger.info('No model_2d_index attribute found, creating index ...')
            self.get_model_2d_index()

        cellidx = self.model_2d_index(xy)
        if indices is None:
            coupled_indices = {i: idx for i, idx in enumerate(cellidx)}
        else:
            coupled_indices = {i: idx for i, idx in zip(indices, cellidx)}
        return coupled_indices, dictinvert(coupled_indices)


class PCR_model(_model):
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
        super(PCR_model, self).__init__(pcr_bmi, config_fn, 'PCRGLOB-WB', 'day',
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

    def get_model_2d_index(self):
        """Get PCR model bounding box, resolution and index based on the
        landmask map file.

        This function creates the model_2d_index attribute function
        which can be used to find the coresponding cell to x, y coordinates.
        """

        fn_map = join(self.model_config['globalOptions']['inputDir'],
                      self.model_config['globalOptions']['landmask'])
        self._landmask_fn = fn_map
        if not isfile(fn_map):
            raise IOError('landmask file not found')
        with rasterio.open(fn_map, 'r') as ds:
            self._model_index = ds.index
            self.model_grid_res = ds.res
            self.model_grid_bounds = ds.bounds
            self.model_grid_shape = ds.shape
            # function for grid row col index
            def model_2d_index(xy, **kwargs):
                r, c = ds.index(*zip(*xy), **kwargs)
                r = np.array(r).astype(int)
                c = np.array(c).astype(int)
                return zip(r, c)
        self.model_2d_index = model_2d_index

        # NOTE: for now decided to keep tuples with index instead linear index.
        # To not have to tranform between both indices all the time.
        #
        # go from tuple with r, c list to linear index. can be added to model_2d_index
        # np.ravel_multi_index((r, c), dims=ds.shape)
        #
        # go from linear index to tuple with r, c lists
        #     def model_index_2_rc(indices):
        #         return np.unravel_index(indices, dims=ds.shape)
        # self.model_index_2_rc = model_index_2_rc

    def get_delta_water(self):
        """get total water volume (discharge, runoff and top water layer) per
        cell per timestep dt"""
        #- retrieve data from PCR-GLOBWB
        current_discharge_pcr  = self.get_var('discharge')
        current_runoff_pcr     = self.get_var('landSurfaceRunoff')
        current_waterlayer_pcr = self.get_var('topWaterLayer')

        # average discharge flux [m3/s]
        water_volume_PCR_rivers = current_discharge_pcr * self.options['tscale'] #sec/dt
        # runoff state [m]
        water_volume_PCR_runoff = current_runoff_pcr * self.get_var('cellArea')
        # topwaterlayer state [m]
        water_volume_PCR_waterlayer = current_waterlayer_pcr * self.get_var('cellArea')
        # sum volumes [m3]
        total_water_volume = water_volume_PCR_rivers + water_volume_PCR_runoff + water_volume_PCR_waterlayer

        return total_water_volume

    def deactivate_LDD(self, index):
        """Deactive LDD at indices

        Parameters
        ----------
        index : list of tuples, str
          list with (x, y) tuples or 'all' to deactivate total grid
        """
        landmask_mv = self.options['landmask_mv']
        # retrieving current LDD map
        if index != 'all':
            # replace LDD values within the hydrodynamic model domain to 5 (pit cell)
            n = len(index) if isinstance(index, list) else len(index[0]) # else: asume list of tuples
            LDD_PCR_new = np.ones(n, dtype=int) * 5
            self.set_var_index(('routing', 'lddMap'), index, LDD_PCR_new, parse_missings=False)
        else: # replace all
            LDD_PCR_new = np.copy(self.get_var(('routing', 'lddMap'), parse_missings=False))
            LDD_PCR_new = np.where(LDD_PCR_new != landmask_mv, 5, landmask_mv)
            # overwriting current with new LDD information
            self.set_var(('routing', 'lddMap'), LDD_PCR_new, parse_missings=False)

class CMF_model(_model):
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
        super(CMF_model, self).__init__(cmf_bmi, config_fn, 'CaMa-Flood', 'sec',
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
                        shutil.copy(
                                src_fn, dst_path)
                    self.update_config(sec, opt, '"./{}/{}"'.format(folder, fn))

    def set_inpmat_file(self, bounds, res, olat='NtoS'):
        """Set the CMF inpmat file model based on the grid definition of upstream
        model"""
        if not abs(res[0]) == abs(res[1]):
            raise ValueError('lat and lon resolution should be the same in regular grid')
        westin  = bounds.left
        eastin  = bounds.right
        northin = bounds.top
        southin = bounds.bottom
        # generate inpmat
        ddir = self.model_data_dir
        msg2 = './generate_inpmat {} {} {} {} {} {:s}'.format(
                        abs(res[0]), westin, eastin, northin, southin, olat)
        logger.info(msg2)
        subcall(msg2, cwd=ddir)
        # set new inpmat and diminfo in config
        rel_path = os.path.relpath(dirname(self.config_fn), ddir)
        inpmatOptions = {'NINPUT': {'CINPMAT': '"{:s}/inpmat-tmp.bin"'.format(rel_path),
                                     'LBMIROF': '.TRUE.'
                                    },
                        'NMAP': {'CDIMINFO': '"{:s}/diminfo_tmp.txt"'.format(rel_path)
                                },
                        'NCONF': {'DROFUNIT': '1'} #  SI units [m]
                        }
        self.set_config(inpmatOptions)

    def get_model_2d_index(self, fn=None):
        # get input file
        if fn is None:
            path = join(self.model_data_dir, 'hires', '*.catmxy.tif')
            fns = glob.glob(path)
            if len(fns) == 0:
                raise IOError("catmxy.tif file not found")
            elif len(fns) > 1:
                raise NotImplemented("Not yet implemented for mulitple regions")
            fn = fns[0]
        else:
            if not isfile(fn):
                raise IOError("catmxy.tif file not found")
        # set model grid paramters
        with rasterio.open(fn, 'r') as ds:
            ncount = ds.count
            assert ncount == 2, "{:s} file should have two layers".format(fn)
            self._catmxy_fn = fn
            self.model_grid_res = ds.res
            self.model_grid_bounds = ds.bounds
            self.model_grid_shape = ds.shape

        # define index function.
        # read catmxy temporary into memory when this function is called
        def model_2d_index(xy, **kwargs):
            with rasterio.open(fn, 'r', driver='GTiff') as ds:
                nrows, ncols = ds.shape
                rows, cols = ds.index(*zip(*xy))
                rows, cols = np.atleast_1d(rows).astype(int), np.atleast_1d(cols).astype(int)
                invalid = np.logical_and.reduce((rows<0,rows>=nrows,cols<0,cols>=ncols))
                if np.any(invalid):
                    raise IndexError('XY coordinates outside of CMF domain')
                cmf_idx = ds.read()[:, rows[:, None], cols[:, None]].squeeze()
                cmf_cols, cmf_rows = cmf_idx
            return zip(cmf_rows, cmf_cols)
        self.model_2d_index = model_2d_index

    def get_var(self, name, parse_missings=True, *args, **kwargs):
        var = super(CMF_model, self).get_var(name, parse_missings=parse_missings)
        # return var with switched axis (fortran to python translation)
        return var.reshape(var.shape[::-1])

    def get_current_time(self):
        t = super(CMF_model, self).get_current_time()
        return CMFtime_2_datetime(t)

    def get_start_time(self):
        t = super(CMF_model, self).get_start_time()
        return CMFtime_2_datetime(t)


class DFM_model(_model):
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
        super(DFM_model, self).__init__(dfm_bmi, config_fn, 'Delft3D-FM', 'sec',
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

    def get_model_coords(self):
        """Get DFM model coordinates for 1D and 2D mesh via BMI. The DFM model
        should be initialized first in order to access the variables."""

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
        # 1d coords
        n1d = len(self.model_1d_coords)
        self.model_1d_indices = np.arange(n1d, dtype=np.int32) + self._1d2d_idx
        # build spatial rtree index of points2
        self.model_1d_rtree = rtree.index.Index()
        for i, xy in enumerate(self.model_1d_coords):
            self.model_1d_rtree.insert(i+n1d, xy) # return index including 2d
        def model_1d_index(xy, n=1):
            if isinstance(xy, tuple):
                xy = [xy]
            return [list(self.model_1d_rtree.nearest(xy0, 1))[0] for xy0 in xy]
        self.model_1d_index = model_1d_index

    def get_model_2d_index(self):
        """Creat a spatial index for the 2d mesh center coordinates.
        A model_2d_index attribute funtion is created to find the nearest
        2d cell center"""
        if not hasattr(self, 'model_2d_coords'):
            self.get_model_coords()

        # 2d coords
        if not only_1d:
            # build spatial rtree index of points2
            self.model_2d_rtree = rtree.index.Index()
            for i, xy in enumerate(self.model_2d_coords):
                self.model_2d_rtree.insert(i, xy)
            def model_2d_index(xy, n=1):
                if isinstance(xy, tuple):
                    xy = [xy]
                return [list(self.model_2d_rtree.nearest(xy0, 1))[0] for xy0 in xy]
            self.model_2d_index = model_2d_index

    # def calculateDeltaWater(self, total_water_volume, af_list):
    #     delta_water = []
    #     for i in xrange(len(total_water_volume)):
    #         added_water_level = total_water_volume[i] / af_list[i]
    #         delta_water.extend(added_water_level)


# utils
def dictinvert(d):
    inv = {}
    for k, v in d.iteritems():
        keys = inv.setdefault(v, [])
        keys.append(k)
    return inv

def CMFtime_2_datetime(t):
    # internal CMF time is in minutes since 1850
    return datetime(1850, 1, 1) + timedelta(minutes = t)

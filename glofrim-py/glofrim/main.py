# -*- coding: utf-8 -*-
# TODO: write main docstring
# TODO: add decorator to child functions to get docstring from parent function

import re
import logging
import warnings
import shutil
import string
import os
from os import mkdir
from os.path import isdir, join, basename, dirname, abspath, isfile, isabs
import numpy as np
from bmi.wrapper import BMIWrapper

# local libraries
from utils import config_to_dict, dict_to_config, ConfigParser

log_fmt = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt, filemode='w')
logger = logging.getLogger(__name__)

# wrapper around BMI
class BMI_model_wrapper(object):
    def __init__(self, bmi, config_fn, name, t_unit,
                 model_data_dir, forcing_data_dir, out_dir,
                 options, configparser=ConfigParser(), si_unit_conversions={}, **kwargs):
        """initialize the model BMI class and model configuration file"""
        # TODO: extend docstring
        self.bmi = bmi
        self.name = name # string used
        self.t_unit = t_unit
        # set model paths
        self.config_fn = abspath(config_fn)
        self._configparser = configparser
        self.model_data_dir = abspath(model_data_dir)
        self.forcing_data_dir = abspath(forcing_data_dir)
        self.out_dir = abspath(out_dir)
        if not isdir(self.out_dir):
            mkdir(self.out_dir)
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


    ## model configuration
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
        #addded if-switch for LFP; 31-May-2018; JMH
        if (self.name not in ['LFP']):
            self.model_config = config_to_dict(self.config_fn,
                                           cf=self._configparser,
                                           **kwargs)
                                           
        if (self.name in ['LFP']):
			#- adapted from
			#- https://stackoverflow.com/questions/2819696/parsing-properties-file-in-python/25493615#25493615

			fo = self.config_fn
			f = open(fo, 'rw') # open LFP par-file
			fake_config = '[dummysection]\n' + f.read() # add dummy section header
			print fake_config
										   
			# EXPERIMENTAL #
			tmpFile = os.path.join(os.path.dirname(self.config_fn), 'tmp.par') # create tmp-file
			file = open(tmpFile, 'w+')
			file.write(fake_config)
			file.close()
			
			self.config_fn = tmpFile
			
			#- PROBLEM: parsing only works with "=" signs...
			
			self.model_config = config_to_dict(self.config_fn,
                                           cf=self._configparser,
                                           **kwargs) # create dictionary
										   
			#os.remove(tmpFile)
                                           

    def write_config(self, **kwargs):
        """The internal model_config dictionary is written to the out_dir. This
        step should be excecuted just before the model initialization."""
        self.config_fn = join(self.out_dir, basename(self.config_fn))
        #addded if-switch for LFP; 31-May-2018; JMH
        if (self.name not in ['LFP']):
			dict_to_config(self.model_config, self.config_fn,
                       cf=self._configparser, **kwargs)
			logger.info('Ini file for {:s} written to {:s}'.format(self.name, self.config_fn))
        if (self.name in ['LFP']):
			print 'LFP section of write_config'
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

    ## model initialization/finalize/spinup
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
        self.bmi.spinup()

    ## exchange states
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

    ## run timestep dt in model
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

    ## time attributes
    def get_start_time(self):
        "get model start (initialization) time"
        return self.bmi.get_start_time()

    def get_current_time(self):
        "get model current time"
        return self.bmi.get_current_time()

    def get_time_step(self):
        "get model current time step"
        return self.bmi.get_time_step()

    ## var info
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
        - WFL : regular grid (row, col)
        - CMF : irregular unit catchment grid (row, col)
        - DFM : flexible mesh (flat index)
        - LFP : regular grid (flat index)

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
        if (other.name not in ['DFM', 'LFP']) or (self.name not in ['PCR', 'CMF', 'WFL']):
            msg = 'Grid to 1D coupling has only been implemented for PCR to' + \
                  ' CMF (upstream) to DFM (downstream) coupling'
            raise NotImplementedError(msg)

        if not other.initialized:
            msg = 'The hydrodynamic model should be initialized first to obtain ' + \
                      ' the 1D model coordinates via BMI'
            raise AssertionError(msg)
        if not hasattr(other, 'model_1d_coords'):
            other.get_model_coords()
        area_other = other.get_area_1d()

        logger.info('Coupling {} grid to {} 1D nodes.'.format(self.name, other.name))
        # get cell indices at 1D coordinates
        cellidx, valid = self.model_2d_index(other.model_1d_coords)
        rows, cols = zip(*cellidx)
        rows, cols = np.atleast_1d(rows), np.atleast_1d(cols)
        # set indices to easily exchange variables
        if not np.all(valid):
            logger.warning('1D nodes found outside of valid 2D domain')
        other.coupled_idx = other.model_1d_indices[valid]
        self.coupled_idx = (rows[valid], cols[valid]) # tuple of (row, col) arrays
        # create coupled 1-to-1 downstream-to-upstream indices dictionary
        other.coupled_dict = {i: (r, c) for i, r, c in zip(other.coupled_idx, *self.coupled_idx)}
        # invert dictionary for 1-to-n upstream-to-downstream coupling
        self.coupled_dict = dictinvert(other.coupled_dict)

        logger.info('Getting fraction of coupled 1d nodes based on area.')
        # get area fractions of coupled cells
        area_frac = {} # [m2/m2]
        for _, idx in list(self.coupled_dict.items()):
            area_other_sum_idx = np.sum(area_other[idx])
            afs = {i: area_other[i]/area_other_sum_idx for i in idx}
            area_frac.update(afs)
        self.coupled_area_frac = np.array([area_frac[i] for i in other.coupled_idx])

        # get mask with 1) coupled runoff and 2) coupled discharge
        rows, cols = zip(*self.coupled_dict.keys())
        rows, cols = np.atleast_1d(rows), np.atleast_1d(cols)
        rows_us, cols_us = self.dd.next_upstream(rows, cols) # find upstream cells where to couple discharge
        # mask of the coupled domain
        coupled_mask = np.zeros(self.model_grid_shape)
        coupled_mask[rows_us, cols_us] = 2 # couple discharge to next downstream cell
        coupled_mask[rows, cols] = 1 # couple runoff
        self.coupled_mask = coupled_mask
        pass

    def couple_grid_to_grid(self, other):
        """Couple external grid to internal model 2d grid.

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

        # deactivate routing in full domain
        self.deactivate_routing('all')

# utils
def dictinvert(d):
    inv = {}
    for k, v in d.iteritems():
        keys = inv.setdefault(v, [])
        keys.append(k)
    return inv

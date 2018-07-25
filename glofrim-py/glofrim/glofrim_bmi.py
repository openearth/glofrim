from collections import OrderedDict
import os
from os.path import join, isfile, abspath, dirname, basename, normpath
from datetime import datetime, timedelta
from configparser import ConfigParser
import logging
import numpy as np

from utils import setlogger
from gbmi import EBmi
import glofrim_lib as glib 
from pcr_bmi import PCR
from cmf_bmi import CMF
from dfm_bmi import DFM
from wfl_bmi import WFL
from lfp_bmi import LFP
from spatial_coupling import SpatialCoupling, groupby_sum

class Glofrim(EBmi):

    """
    csdms BMI implementation runner for combined hydraulic and hydrological models
    """
    _name = 'GLOFRIM'
    _version = '2.0'
    _models = {'PCR': PCR, 'CMF': CMF, 'DFM': DFM, 'WFL': WFL, 'LFP': LFP}
    _init_pre_exchange = ['DFM'] # need to initialize before we can know the grid

    def __init__(self):
        self.bmimodels = OrderedDict()
        self.exchanges = []
        self._var_sep = "."
        self._mult_sep = "*"
        self._ind_sep = "@"

        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False
        self.obs = None

    def _check_long_var_name(self, long_var_name):
        if not (len(long_var_name.split(self._var_sep)) == 2):
            msg = 'use "model{}var" syntax for long_var_name'.format(self._var_sep)
            raise ValueError(msg)
        return long_var_name.split(self._var_sep)

    def _check_initialized(self):
        self.initialized = np.all([self.bmimodels[mod].initialized for mod in self.bmimodels])
        return self.initialized

    def _check_dt(self):
        for mod in self.bmimodels:
            dt_mod = self.bmimodels[mod].get_time_step()
            dt_mod = (dt_mod.days*86400 + dt_mod.seconds)
            _dt = (self._dt.days*86400 + self._dt.seconds)
            dt_diff = _dt % dt_mod
            if dt_diff > 0:
                msg = 'A whole number of {:s} model timesteps ({:d} sec) should fit in the GLOFRIM exchange time step ({:d} sec). ' \
                      'Please change the model timestep in the model configuration file'  
                self.logger.error(msg)
                raise AssertionError(msg.format(mod, dt_mod, _dt))
        return self.initialized

    """
    Model Control Functions
    """
    def initialize_config(self, config_fn, env_fn=None):       
        if self.initialized:
            raise Warning("model already initialized, it's therefore not longer possible to initialize the config")
        # config settings
        self._config_fn = abspath(config_fn)
        self._root = dirname(config_fn)
        self._config = glib.configread(self._config_fn, encoding='utf-8', 
                                       cf=ConfigParser(inline_comment_prefixes=('#')))
       
        # environment file -> merge with config if given
        if env_fn is not None and isfile(abspath(env_fn)):
            env_fn = abspath(env_fn)
            env = glib.configread(env_fn, encoding='utf-8', 
                                  cf=ConfigParser(inline_comment_prefixes=('#')))
            for sect in env.sections():
                if sect not in self._config.sections():
                    self._config.add_section(sect)
                for opt in env.options(sect):
                    self._config.set(sect, opt, glib.getabspath(env.get(sect, opt), dirname(env_fn)))
        if self._config.has_option('models', 'root_dir'):
            self._root = glib.getabspath(self._config.get('models', 'root_dir'), self._root)
            self._config.remove_option('models', 'root_dir')
        
        ## parse glofrim config
        # initialize bmi component and it's configuration 
        if not self._config.has_section('models'):
            raise ValueError('GLOFRIM ini misses a "models" section')
        for mod in self._config.options('models'):
            bmi_kwargs = {}
            _bmi = self._models[mod]
            # check if bmi component requires engine
            if 'engine' in _bmi.__init__.__code__.co_varnames:
                if not self._config.has_option('engines', mod):
                    msg = 'GLOFRIM ini or environment file misses a "engines" section with {} option'
                    raise ValueError(msg.format(mod))
                engine_path = glib.getabspath(self._config.get('engines', mod), self._root)
                bmi_kwargs = dict(engine = engine_path)
            # initialize bmi component
            self.bmimodels[mod] = _bmi(**bmi_kwargs)
            # initialize config of bmi component
            modconf = glib.getabspath(self._config.get('models', mod), self._root)
            self.bmimodels[mod].initialize_config(modconf)

        # parse echanges section
        self.exchanges = self.set_exchanges()

        # combined model time
        self._dt = timedelta(seconds=int(self._config.get('coupling' ,'dt', fallback=86400)))
        # set start_time & end_time if given
        if self._config.has_option('coupling', 'start_time'):
            start_time = datetime.strptime(self._config.get('coupling', 'start_time'), "%Y-%m-%d")
            self.set_start_time(start_time)
        if self._config.has_option('coupling', 'end_time'):
            end_time = datetime.strptime(self._config.get('coupling', 'end_time'), "%Y-%m-%d")
            self.set_end_time(end_time)
        self._timeunit = 'sec'
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        self.initialized=False
        # check model dt
        self._check_dt()

    def set_exchanges(self):
        # parse exchanges
        if not self._config.has_section('exchanges'):
            raise ValueError('GLOFRIM ini misses a "exchanges" section')
        self.exchanges = []
        model_called, vars_set = [], []
        for ex_from in self._config.options('exchanges'):
            ex_to = self._config.get('exchanges', ex_from)
            # break up starting from back
            exdict = {}
            ## FROM
            # model
            if self._ind_sep in ex_from:
                ex_from, ind_from = ex_from.split(self._ind_sep)
            else:
                ind_from = 'grid'
            ex_from = ex_from.split(self._var_sep)
            mod_from, vars_from = ex_from[0], '.'.join(ex_from[1:])
            if not (mod_from in self.bmimodels.keys()):
                raise ValueError("unkown model name {}".format(mod_from))
            exdict['from_mod'] = mod_from
            from_bmi =  self.bmimodels[mod_from]
            # variables
            exdict['from_vars'], exdict['from_unit'] = [], []
            vars_from = vars_from.split('*')
            from_model_vars = from_bmi._input_var_names + from_bmi._output_var_names
            for i, var_from in enumerate(vars_from):
                # last var may be a scalar multiplier
                if i > 0 and i == (len(vars_from) - 1):
                    try:
                        var_from = float(var_from)
                    except ValueError:
                        pass
                exdict['from_vars'].append(var_from) # model variable
                exdict['from_unit'].append(from_bmi.get_var_units(var_from))
                if isinstance(var_from, str) and var_from not in from_model_vars:
                    self.logger.warning("Unkonwn variable {} for model {}".format(var_from, from_bmi._name))
            ## TO
            # index
            if self._ind_sep in ex_to:
                ex_to, ind_to = ex_to.split(self._ind_sep)
            else:
                ind_to = 'grid'
            # model
            ex_to = ex_to.split(self._var_sep)
            mod_to, vars_to = ex_to[0], '.'.join(ex_to[1:])
            if not (mod_to in self.bmimodels.keys()):
                raise ValueError("unkown model name {}".format(mod_to))
            to_bmi = self.bmimodels[mod_to]
            exdict['to_mod'] = mod_to
            # variables
            exdict['to_vars'],  exdict['to_unit'] = [], []
            vars_to = vars_to.split('*')
            to_model_vars = to_bmi._input_var_names + to_bmi._output_var_names
            for i, var_to in enumerate(vars_to):
                # last var may be a scalar multiplier
                if i > 0 and i == (len(vars_to) - 1):
                    try:
                        var_to = float(var_to)
                    except ValueError:
                        pass
                exdict['to_vars'].append(var_to) # model variable
                exdict['to_unit'].append(to_bmi.get_var_units(var_to))
                if isinstance(var_to, str) and var_to not in to_model_vars:
                    self.logger.Warning("Unkonwn variable {} for model {}".format(var_to, to_bmi._name))
            # with second call to set variable add instead of replace
            if exdict['to_vars'][0] in vars_set:
                exdict['add'] = True 
            else:
                vars_set.append(exdict['to_vars'][0])
            ## SpatialCoupling
            # create SpatialCoupling object. Actual coupling happens later
            sc_kwargs = dict(to_bmi=to_bmi, from_bmi=from_bmi)
            if isfile(ind_to):
                sc_kwargs.update(filename=ind_to, method='from_file')
            else:
                coupling_method = '{}_2_{}'.format(ind_from, ind_to)
                sc_kwargs.update(method=coupling_method)
            exdict['coupling'] = SpatialCoupling(**sc_kwargs)
            ## UPDATE
            #  update model before first time it is called in mod_from
            if not mod_from in model_called:
                self.exchanges.append(('update', mod_from))
                model_called.append(mod_from)
            self.exchanges.append(('exchange', exdict))
        # add last model updates
        for mod in self.bmimodels.keys():
            if not mod in model_called:
                self.exchanges.append(('update', mod))
                model_called.append(mod)
        return self.exchanges

    def _set_spatial_coupling(self):
        for item in self.exchanges:
            if item[0] == 'exchange':
                item[1]['coupling'].couple()

    def initialize_model(self, **kwargs):
        if not hasattr(self, '_config_fn'):
            raise Warning('run initialize_config before initialize_model')
        # some models (e.g. DFM) need to be intialized before we access it spatial information
        for mod in self.bmimodels:
            if mod in self._init_pre_exchange:
                self.bmimodels[mod].initialize_model()
        # set spatial coupling
        self._set_spatial_coupling()
        # initialize other models (CMF needs to be initialized after spatial coupling!)
        for mod in self.bmimodels:
            if mod not in self._init_pre_exchange:
                self.bmimodels[mod].initialize_model()
        self._check_initialized()
        # sync start and endtimes
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # check model dt
        self._check_dt()
        # TODO set observation points

    def initialize(self, config_fn):
        self.initialize_config(config_fn)
        self.initialize_model()

    def update(self):
        if not self.initialized:
            raise Warning("models should be initialized first")
        if self._t >= self.get_end_time():
		    raise Exception("endTime already reached, model not updated")
        next_t = self.get_current_time() + self.get_time_step()
        for item in self.exchanges:
            if item[0] == 'update':
                self.bmimodels[item[1]].update_until(next_t)
            elif item[0] == 'exchange':
                self.exchange(**item[1])
        self._t = next_t

    def update_until(self, t):
        if (t<self._t) or t>self.get_end_time():
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update()

    def spinup(self):
        """PCR specific spinup function"""
        raise NotImplementedError

    def finalize(self):
        for mod in self.bmimodels:
            self.bmimodels[mod].finalize()

    def exchange(self, from_mod, to_mod, from_vars, to_vars, coupling, add=False, **kwargs):
        self.logger.debug("get data from {}.{}; set to {}.{}".format(from_mod, from_vars[0], to_mod, to_vars[0]))
        if coupling.at_indices():
            self.exchange_at_indices(from_mod, to_mod, from_vars, to_vars, coupling, add=add, **kwargs)
        else:
            self.exchange_same_grid(from_mod, to_mod, from_vars, to_vars, add=add, **kwargs)

    def exchange_same_grid(self, from_mod, to_mod, from_vars, to_vars, add=False, **kwargs):
        # get FROM data & translate
        shape = self.bmimodels[from_mod].get_var_shape(from_vars[0])
        vals = np.ones(shape)
        for from_var in from_vars:
            if isinstance(from_var, float): # scalar multiplier
                vals *= from_var
            else:
                vals *= self.bmimodels[from_mod].get_value(from_var)
        self.logger.debug('total get {:3f}'.format(np.nansum(vals)))
        # get TO data & translate
        for var in to_vars[1:]:
            if isinstance(var, float):
                div = var
            else:
                div = self.bmimodels[to_mod].get_value(var) 
            assert np.all(np.not_equal(div, 0)) # make sure not to divide by zero
            vals /= div
        # SET data
        if add:  # add to current state
            vals += self.bmimodels[to_mod].set_value(to_vars[0], vals) # TO   
        self.bmimodels[to_mod].set_value(to_vars[0], vals)

    def exchange_at_indices(self, from_mod, to_mod, from_vars, to_vars, coupling, add=False, **kwargs):
        """
        exchanges values between variables of two models
        """
        # get FROM data & translate
        vals_get = np.ones(coupling.from_ind.size)
        for from_var in from_vars:
            if isinstance(from_var, float): # scalar multiplier
                vals_get *= from_var
            else:
                vals_get *= self.bmimodels[from_mod].get_value_at_indices(from_var, coupling.from_ind)
        # TRANSFORM data from coupling.from_ind.size > coupling.to_ind.size
        if coupling.from_grp_n.max() > 1: # aggregate data from
            vals_get = groupby_sum(vals_get, coupling.from_grp)
        if coupling.to_grp_n.max() > 1: # split using fractions
            vals_set = np.repeat(vals_get, coupling.to_grp_n) * coupling.get_frac()
        else:
            vals_set = vals_get
        vol_diff = np.nansum(vals_get) - np.nansum(vals_set)
        if vol_diff > 1E-2:
            self.logger.warming('Large difference in water volume from and to: {:.4f} m3'.format(vol_diff))
        self.logger.debug('total get {:3f}'.format(np.nansum(vals_get)))
        # get TO data & translate
        for var in to_vars[1:]:
            if isinstance(var, float):
                div = var
            else:
                div = self.bmimodels[to_mod].get_value_at_indices(var, coupling.to_ind)
            assert np.all(np.not_equal(div, 0)) # make sure not to divide by zero
            vals_set /= div
        # SET data
        if add: # add to current state
            vals_set += self.bmimodels[to_mod].get_value_at_indices(to_vars[0], coupling.to_ind)
        self.bmimodels[to_mod].set_value_at_indices(to_vars[0], coupling.to_ind, vals_set)

    """
    Model Information Functions
    """
    def get_model_type(self):
        return {mod: self.bmimodels[mod].get_model_type() for mod in self.bmimodels}

    def get_component_name(self):
        return {mod: self.bmimodels[mod].get_component_name() for mod in self.bmimodels}

    def get_input_var_names(self):
        return {mod: self.bmimodels[mod].get_input_var_names() for mod in self.bmimodels}

    def get_output_var_names(self):
        return {mod: self.bmimodels[mod].get_output_var_names() for mod in self.bmimodels}


    """
    Variable Information Functions
    """

    def get_var_type(self, long_var_name):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_type(var)

    def get_var_units(self, long_var_name):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_units(var)

    def get_var_rank(self, long_var_name):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_rank(var)

    def get_var_size(self, long_var_name):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_size(var)

    def get_var_shape(self, long_var_name):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_shape(var)

    def get_var_nbytes(self, long_var_name):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_nbytes(var)

    
    def get_start_time(self):
        self._startTime = max([self.bmimodels[mod].get_start_time() for mod in self.bmimodels])
        return self._startTime

    def get_current_time(self):
        return self._t

    def get_end_time(self):
        self._endTime = min([self.bmimodels[mod].get_end_time() for mod in self.bmimodels])
        return self._endTime

    def get_time_step(self):
        
        return self._dt
        
    def get_time_units(self):
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, **kwargs):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_value(var, **kwargs)

    def get_value_at_indices(self, long_var_name, inds, **kwargs):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_value_at_indices(var, inds, **kwargs)

    def set_value(self, long_var_name, src, **kwargs):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].set_value(var, src, **kwargs)

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].set_value_at_indices(var, inds, src, **kwargs)



    """
    Grid Information Functions
    """
    
    def get_grid_type(self):
        return {mod: self.bmimodels[mod].get_grid_type() for mod in self.bmimodels}

    """
    Observation points
    """
    def get_obs(self):
        pass

    def set_obs(self):
        pass

    def index(self, x, y, mod, in1d=False):
        if self.bmimodels[mod].grid is None:
            raise AssertionError("{}: model grid not set".format(mod))
        if in1d:
            if self.bmimodels[mod].grid._1d is None:
                raise AssertionError("{}: model 1d network not set".format(mod))
            idx = self.bmimodels[mod].grid._1d.index(x,y)
        else:
            idx = self.bmimodels[mod].grid.index(x,y)
        return idx
    
    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        for mod in self.bmimodels:
            self.bmimodels[mod].set_start_time(start_time)
        self.get_start_time() # sync start time
        self._t = self._startTime
       
    def set_end_time(self, end_time):
        for mod in self.bmimodels:
            self.bmimodels[mod].set_end_time(end_time)
        self.get_end_time() # sync start time

    def get_attribute_names(self):
        """
        return a list with all config option names per model
        """
        return {mod: self.bmimodels[mod].get_attribute_names() for mod in self.bmimodels}

    def set_out_dir(self, out_dir):
        """
        set and create outdir/model_name/ paths for output file
        """
        for mod in self.bmimodels:
            path = join(out_dir, mod)
            if not os.path.isdir(path):
                os.makedirs(path)
            self.bmimodels[mod].set_out_dir(path)

    def get_attribute_value(self, attribute_name):
        """
        sets attribute value in in underlying model
            attribute_name should use the following convention model_name.section_name:attribute_name
        """
        mod, attr = self._check_long_var_name(attribute_name)
        return self.bmimodels[mod].get_attribute_value(attr)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        """
        sets attribute value in in underlying model
            attribute_name should use the following convention model_name.section_name:attribute_name
        """
        mod, attr = self._check_long_var_name(attribute_name)
        return self.bmimodels[mod].set_attribute_value(attr, attribute_value)

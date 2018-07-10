from collections import OrderedDict
import os
from os.path import join, isfile, abspath, dirname, basename, normpath
from datetime import datetime, timedelta
from configparser import ConfigParser
import logging
import numpy as np

from utils import setlogger
from gbmi import EBmi
import gbmi_lib as glib 
from pcr_bmi import PCR
from cmf_bmi import CMF
from spatial_coupling import SpatialCoupling

class Glofrim(EBmi):

    """
    csdms BMI implementation runner for combined hydraulic and hydrological models
    """
    _name = 'GLOFRIM'
    _version = '2.0'
    _models = {'PCR': PCR, 'CMF': CMF}

    def __init__(self):
        self.bmimodels = OrderedDict()
        self.exchanges = []
        self._var_sep = "."
        self._mult_sep = "*"
        self._ind_sep = "@"

        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False

    def _check_long_var_name(self, long_var_name):
        if not (len(long_var_name.split(self._var_sep)) == 2):
            msg = 'use "model{}var" syntax for long_var_name'.format(self._var_sep)
            raise ValueError(msg)
        return long_var_name.split(self._var_sep)

    def _check_initialized(self):
        self.initialized = np.all([self.bmimodels[mod].initialized for mod in self.bmimodels])
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
                self._config.add_section(sect)
                for opt in env.options(sect):
                    self._config.set(sect, opt, glib.getabspath(env.get(sect, opt), dirname(env_fn)))
        
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

    def set_exchanges(self):
        # parse exchanges
        if not self._config.has_section('exchanges'):
            raise ValueError('GLOFRIM ini misses a "exchanges" section')
        model_called = []
        for ex_from in self._config.options('exchanges'):
            ex_to = self._config.get('exchanges', ex_from)
            if not ((len(ex_from.split(self._var_sep)) == 2) and 
                    (len(ex_to.split(self._var_sep)) == 2)):
                msg = 'use "From_model{}var = To_model{}var" syntax for exchange options'.format(self._var_sep, self._var_sep)
                raise ValueError(msg)
            # break up starting from back
            exdict = {}
            # FROM
            if self._ind_sep in ex_from:
                ex_from, exdict['from_ind'] = ex_from.split(self._ind_sep)
            if self._mult_sep in ex_from:
                ex_from, m = ex_from.split(self._mult_sep)
                exdict['multiplier'] = float(m)
            mod_from, exdict['from_var'] = ex_from.split(self._var_sep)
            if not (mod_from in self.bmimodels.keys()):
                raise ValueError("unkown model name {}".format(mod_from))
            exdict['from_mod'] = mod_from
            # TO
            if self._ind_sep in ex_to:
                ex_to, coupling_method = ex_to.split(self._ind_sep)
            mod_to, exdict['to_var'] = ex_to.split(self._var_sep)
            if not (mod_to in self.bmimodels.keys()):
                raise ValueError("unkown model name {}".format(mod_to))
            exdict['to_mod'] = mod_to
            # do spatial coupling
            exdict['coupling'] = SpatialCoupling()
            if isfile(coupling_method):
                json_fn = coupling_method
                exdict['coupling'].from_file(json_fn)
            elif hasattr(exdict['coupling'], coupling_method):
                getattr(exdict['coupling'], coupling_method)(self.bmimodels[mod_to], self.bmimodels[mod_from])
            else:
                raise AttributeError('Coupling method {} does not exist'.format(exdict['coupling']))
            # update model before first time it is called in mod_from
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

    def initialize_model(self, **kwargs):
        self.set_exchanges()
        if not hasattr(self, '_config_fn'):
            raise Warning('run initialize_config before initialize_model')
        for mod in self.bmimodels:
            self.bmimodels[mod].initialize_model()
        self._check_initialized()

    def initialize(self, config_fn):
        self.initialize_config(config_fn)
        self.initialize_model()

    def update(self):
        if not self.initialized:
            raise Warning("models should be initialized first")
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        next_t = self.get_current_time() + self.get_time_step()
        for item in self.exchanges:
            if item[0] == 'update':
                self.bmimodels[item[1]].update_until(next_t)
            elif item[0] == 'exchange':
                self.exchange(**item[1])
        self._t = next_t

    def update_until(self, t):
        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update()

    def spinup(self):
        """PCR specific spinup function"""
        raise NotImplementedError

    def finalize(self):
        for mod in self.bmimodels:
            self.bmimodels[mod].finalize()

    def exchange(self, from_mod, to_mod, from_var, to_var, multiplier=1., coupling=SpatialCoupling()):
        """
        exchanges values between variables of two models
        """
        from_ind = coupling.from_ind
        to_ind = coupling.to_ind
        # TODO deal with dict
        if len(from_ind) == 0:
            value = self.bmimodels[from_mod].get_value(from_var) * multiplier # FROM
            self.bmimodels[to_mod].set_value(to_var, value) # TO
        else:
            value = self.bmimodels[from_mod].get_value_at_indices(from_var, from_ind) * multiplier # FROM
            self.bmimodels[to_mod].set_value_at_indices(to_var, value, to_ind) # TO
        
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
        return max([self.bmimodels[mod].get_start_time() for mod in self.bmimodels])
    
    def get_current_time(self):
        return self._t

    def get_end_time(self):
         return min([self.bmimodels[mod].get_end_time() for mod in self.bmimodels])

    def get_time_step(self):
        #NOTE get_time_step in PCR bmi returns timestep as int instead of dt
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
        return self.bmimodels[mod].get_value_at_indices(var, inds, src, **kwargs)

    """
    Grid Information Functions
    """
    
    def get_grid_type(self):
        return {mod: self.bmimodels[mod].get_grid_type() for mod in self.bmimodels}

    def get_grid_transform(self):
        return {mod: self.bmimodels[mod].get_grid_type() for mod in self.bmimodels
                if self.bmimodels[mod].get_grid_type == 2} # only for RECTILINEAR grids

    def get_grid_shape(self):
        return {mod: self.bmimodels[mod].get_grid_shape() for mod in self.bmimodels}

    def get_grid_bounds(self):
        return {mod: self.bmimodels[mod].get_grid_bounds() for mod in self.bmimodels}

    def get_grid_res(self):
        return {mod: self.bmimodels[mod].get_grid_res() for mod in self.bmimodels}

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        for mod in self.bmimodels:
            self.bmimodels[mod].set_start_time(start_time)
       
    def set_end_time(self, end_time):
        for mod in self.bmimodels:
            self.bmimodels[mod].set_end_time(end_time)

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


from collections import OrderedDict
import os
from os.path import join, isfile, abspath, dirname, basename, normpath
from datetime import datetime, timedelta
from configparser import ConfigParser
import logging
import numpy as np

from glofrim.utils import setlogger, closelogger, add_file_handler
from glofrim.gbmi import EBmi
import glofrim.glofrim_lib as glib
from glofrim.pcr_bmi import PCR
from glofrim.cmf_bmi import CMF
from glofrim.dfm_bmi import DFM
from glofrim.wfl_bmi import WFL
from glofrim.lfp_bmi import LFP
from glofrim.sfincs_bmi import Sfincs
from glofrim.spatial_coupling import SpatialCoupling, groupby_sum

class Glofrim(EBmi):
    """Central CSDMS-compliant BMI implementation; 
    generic script providing function sceleton to other model-specified BMIs; 
    defining main and central model and BMI information
    
    Arguments:
        EBmi {class} -- Interface (abstract base class) for a model that implements the CSDMS BMI (Basic Model Interface).
    
    Raises:
        ValueError -- Raised if syntax requesting variable name is not correct
        AssertionError -- Raised if model dt is not a whole number of GLOFRIM dt
        Warning -- Warning is raised if two-step model initialization is not correctly executed
        ValueError -- Raised if config-files of models to be coupled are not defined in ini-file
        ValueError -- Raised if engines (ie executables) for models to be coupled are not specified in ini/env-file
        ValueError -- Raised if no exchanges are specified in GLOFRIM ini-file
        ValueError -- Raised if an unspported model is specified
        ValueError -- Raised if unknown variables are specified
        Warning -- Warning is raised if two-step model initialization is not correctly executed
        Warning -- Raised in models are not initialized before updating
        Exception -- Raised if model end time is already reached; no further updating possible
        Exception -- Raised if specified time is smaller than time step or later than model end time
        NotImplementedError -- SpinUp is not implemented for central BMI
        AssertionError -- Raised if no grid is set before
        AssertionError -- Raised if no 1D network is set before
    """

    _name = 'GLOFRIM'
    _version = '2.0'
    _models = {'PCR': PCR, 'CMF': CMF, 'DFM': DFM, 'WFL': WFL, 'LFP': LFP, 'Sfincs': Sfincs}
    _init_pre_exchange = ['DFM'] # need to initialize before we can know the grid

    def __init__(self, loglevel=logging.INFO):
        self.bmimodels = OrderedDict()
        self.exchanges = []
        self._var_sep = "."
        self._mult_sep = "*"
        self._ind_sep = "@"
        self._coord_sep = "|"

        self._loglevel = loglevel
        self.logger = setlogger(None, self._name, thelevel=loglevel)
        self.wb_logger = setlogger(None, 'wb', thelevel=loglevel, show_in_console=False)
        self.initialized = False
        self.obs = None

    def _check_long_var_name(self, long_var_name):
        #TODO: what's the use of this function? You provide a name and get it back split?
        """Returns variable name of a user-specified variable.
        
        Arguments:
            long_var_name {str} -- variable long name
        
        Raises:
            ValueError -- Raised if syntax requesting variable name is not correct
        
        Returns:
            str -- variable name of specified variable
        """

        if not (len(long_var_name.split(self._var_sep)) == 2):
            msg = 'use "model{}var" syntax for long_var_name'.format(self._var_sep)
            self.logger.error(msg)
            raise ValueError(msg)
        return long_var_name.split(self._var_sep)

    def _check_initialized(self):
        """Checks whether model was already initialized previously.
        criterion needed since some functions require the model to be in initialized state
        
        Returns:
            boolean -- Boolean value whether model is already initialized
        """

        self.initialized = np.all([self.bmimodels[mod].initialized for mod in self.bmimodels])
        return self.initialized

    def _check_dt(self):
        """Checks whether the dt specified for a model is not a whole number of GLOFRIM dt.
        criterion is needed to ensure updates between models remain in phase
        
        Raises:
            AssertionError -- Raised if model dt is not a whole number of GLOFRIM dt
        """

        for mod in self.bmimodels:
            dt_mod = self.bmimodels[mod].get_time_step()
            if not glib.check_dts_divmod(self._dt, dt_mod):
                msg = f"Invalid value for dt {self._dt} in comparison to model dt {dt_mod}. Make sure a whole number of model timesteps fit in the set GLOFRIM dt"
                self.logger.error(msg)
                raise AssertionError(msg)

    ###
    ### Model Control Functions
    ###

    def initialize_config(self, config_fn, env_fn=None):    
        """Initializing the model configuration file.
        aligning GLOFRIM specifications for coupled runs to be consistent with overall run settings; 
        with environment file (env-file) local paths to model engines can be defined
                
        Arguments:
            config_fn {str} -- path to model configuration file
        
        Keyword Arguments:
            env_fn {str} -- path to environment file (default: {None})
        
        Raises:
            Warning -- Warning is raised if two-step model initialization is not correctly executed
            ValueError -- Raised if config-files of models to be coupled are not defined in ini-file
            ValueError -- Raised if engines (ie executables) for models to be coupled are not specified in ini/env-file
        """

        # log to file
        add_file_handler(self.logger, abspath(config_fn).replace('.ini', '.log'))
        if self.initialized:
            msg = "model already initialized, it's therefore not longer possible to initialize the config"
            self.logger.warn(msg)
            raise Warning(msg)
        # config settings
        self._config_fn = abspath(config_fn)
        self._root = dirname(config_fn)
        self.logger.info('Reading ini file..')
        self._config = glib.configread(self._config_fn, encoding='utf-8', 
                                       cf=ConfigParser(inline_comment_prefixes=('#')))
       
        # environment file -> merge with config if given
        if env_fn is not None and isfile(abspath(env_fn)):
            env_fn = abspath(env_fn)
            self.logger.info('Reading env file..')
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
            msg = 'GLOFRIM ini misses a "models" section'
            self.logger.error(msg)
            raise ValueError(msg)
        for mod in self._config.options('models'):
            bmi_kwargs = {'logger': self.logger, 'loglevel': self._loglevel}
            _bmi = self._models[mod]
            # check if bmi component requires engine
            if 'engine' in _bmi.__init__.__code__.co_varnames:
                if not self._config.has_option('engines', mod):
                    msg = 'GLOFRIM ini or environment file misses a "engines" section with {} option'
                    self.logger.error(msg)
                    raise ValueError(msg.format(mod))
                engine_path = glib.getabspath(self._config.get('engines', mod), self._root)
                bmi_kwargs.update(engine = engine_path)
            # read proj string from config file
            if not self._config.has_option('coupling', mod):
                msg = 'GLOFRIM ini or environment file misses a "coupling" section {} option for projection string, assuming lat-lon'
                self.logger.error(msg)
                # raise ValueError(msg.format(mod))
            crs = self._config.get('coupling', mod, fallback='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

            # initialize bmi component
            self.bmimodels[mod] = _bmi(**bmi_kwargs)
            # initialize config of bmi component
            modconf = glib.getabspath(self._config.get('models', mod), self._root)
            self.bmimodels[mod].initialize_config(modconf)
            # initialize grid
            self.bmimodels[mod].get_grid()
            # if grid does not have a crs, it will be taken from the config file
            if self.bmimodels[mod].grid.crs is None:
                self.bmimodels[mod].grid.crs = crs

        # parse exchanges section
        self.exchanges = self.set_exchanges()
        # create logfile for exchange volumes
        add_file_handler(self.wb_logger, abspath(config_fn).replace('.ini', '.wb'), formatter=logging.Formatter("%(message)s"))
        self._wb_header = ['time']
        for mod in self.bmimodels:
            if hasattr(self.bmimodels[mod], '_get_tot_volume_in'):
                self._wb_header.append('{}_tot_in'.format(mod))
            if hasattr(self.bmimodels[mod], '_get_tot_volume_out'):
                self._wb_header.append('{}_tot_out'.format(mod))
        self._wb_header += [item[1]['name'] for item in self.exchanges if item[0] == 'exchange']
        self.wb_logger.info(', '.join(self._wb_header))
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
        """Defines variable exchanges between models as well as unit and time conversions. 
        exchanges can be defined in ini/env-file;
        it is important that all fluxes/states in ini/env-file are to m3 for mass balance.
        
        Raises:
            ValueError -- Raised if no exchanges are specified in GLOFRIM ini-file
            ValueError -- Raised if an unspported model is specified
            ValueError -- Raised if unknown variables are specified
        
        Returns:
            list -- list describing all models, states, and fluxes that will be exchanged
        """

        self.logger.info('Parsing exchanges..')
        # parse exchanges
        if not self._config.has_section('exchanges'):
            msg = 'GLOFRIM ini misses a "exchanges" section'
            self.logger.error(msg)
            raise ValueError(msg)
        self.exchanges = []
        model_called, vars_set = [], []
        to_coords = None
        for ex_i, ex_from in enumerate(self._config.options('exchanges')):
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
                msg = "unkown model name {}".format(mod_from)
                self.logger.error9(msg)
                raise ValueError(msg)
            exdict['from_mod'] = mod_from
            from_bmi =  self.bmimodels[mod_from]
            # variables
            exdict['from_vars'], exdict['from_unit'] = [], []
            vars_from = vars_from.split('*')
            from_model_vars = from_bmi._input_var_names + from_bmi._output_var_names
            for i, var_from in enumerate(vars_from):
                # last var may be a scalar multiplier
                if i > 0: # and i == (len(vars_from) - 1):
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
                # check if manual set coordinates are found
                if self._coord_sep in ind_to:
                    ind_to, to_coords = ind_to.split(self._coord_sep)
                    to_coords = eval(to_coords)
            else:
                ind_to = 'grid'
            # model
            ex_to = ex_to.split(self._var_sep)
            mod_to, vars_to = ex_to[0], '.'.join(ex_to[1:])
            if not (mod_to in self.bmimodels.keys()):
                msg = "unkown model name {}".format(mod_from)
                self.logger.error(msg)
                raise ValueError(msg)
            to_bmi = self.bmimodels[mod_to]
            exdict['to_mod'] = mod_to
            # variables
            exdict['to_vars'],  exdict['to_unit'] = [], []
            vars_to = vars_to.split('*')
            to_model_vars = to_bmi._input_var_names + to_bmi._output_var_names
            for i, var_to in enumerate(vars_to):
                # last var may be a scalar multiplier
                if i > 0: # and i == (len(vars_to) - 1):
                    try:
                        var_to = float(var_to)
                    except ValueError:
                        pass
                exdict['to_vars'].append(var_to) # model variable
                exdict['to_unit'].append(to_bmi.get_var_units(var_to))
                if isinstance(var_to, str) and var_to not in to_model_vars:
                    self.logger.warning("Unkonwn variable {} for model {}".format(var_to, to_bmi._name))
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
                sc_kwargs.update(method=coupling_method, to_coords=to_coords)
            exdict['coupling'] = SpatialCoupling(**sc_kwargs)
            exdict['name'] = '{:s}.{:s}_to_{:s}.{:s}'.format(mod_from, exdict['from_vars'][0], mod_to, exdict['to_vars'][0])
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
        """Couples model grids.
        for one-way coupling, this is a one-to-many coupling, ie one origin cell is coupled
        to multiple destination cells
        """

        for item in self.exchanges:
            if item[0] == 'exchange':
                self.logger.info('set coupled grids for: {}'.format(str(item[1]['coupling'])))
                item[1]['coupling'].couple()

    def initialize_model(self, **kwargs):
        """Initializes the model, ie loading all files and checking for consistency. 
        can only be executed after the model-specific config-file was initialized;
        essential for a successful two-step model initialization and aligned model data.
        
        Raises:
            Warning -- Warning is raised if two-step model initialization is not correctly executed
        """

        if not hasattr(self, '_config_fn'):
            msg = 'run initialize_config before initialize_model'
            self.logger.warn(msg)
            raise Warning(msg)
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
        """Initializes the model following a two-step initialization procedure.
        first, the config-file is initialized and where necessary aligned with overall model settings (e.g. output-dir, start and end time);
        second, the model is actually initialized.
        
        Arguments:
            config_fn {str} -- path to model configuration file
        """

        self.initialize_config(config_fn)
        self.initialize_model()

    def update(self, dt=None):
        """Updating model for a certain time step interval (default: None).
        checks whether model end time is already reached;
        requires that model-specific time step is a whole number of update time step.
        
        Keyword Arguments:
            dt {int} -- update time step interval [sec] at which information is exchanged between models (default: {None})
        
        Raises:
            Warning -- Raised in models are not initialized before updating
            Exception -- Raised if model end time is already reached; no further updating possible
        """

        if not self.initialized:
            msg = "models should be initialized first"
            self.logger.warn(msg)
            raise Warning(msg)
        if self._t >= self.get_end_time():
            msg = "endTime already reached, model not updated"
            self.logger.warn(msg)
            raise Exception(msg)
        # update all models with combined model dt
        wb_dict = {'time': str(self._t)}
        dt = self._dt.total_seconds() if dt is None else dt
        t_next = self._t + timedelta(seconds=dt)
        for item in self.exchanges:
            if item[0] == 'update':
                # LFP deviates from the set timestep using an adaptive timestep if dt is set to large
                # calculate the dt to get to the next timestep
                # NOTE we use "update" instead of "update_until" to getter better logging.
                imod = item[1]
                dt_mod = (t_next - self.bmimodels[imod]._t).total_seconds()
                self.bmimodels[imod].update(dt=dt_mod)
                # get volume totals in and out if the bmi object has this funtion, ortherwise return -9999
                tot_volume_in = getattr(self.bmimodels[imod], '_get_tot_volume_in', lambda: -9999.)()
                wb_dict.update({'{}_tot_in'.format(imod): '{:.2f}'.format(tot_volume_in)})
                tot_volume_out = getattr(self.bmimodels[imod], '_get_tot_volume_out', lambda: -9999.)()
                wb_dict.update({'{}_tot_out'.format(imod): '{:.2f}'.format(tot_volume_out)})
            elif item[0] == 'exchange':
                tot_volume = self.exchange(**item[1])
                wb_dict.update({item[1]['name']: '{:.2f}'.format(tot_volume)})
        # write water balance volumes to file according to header
        self.wb_logger.info(', '.join([wb_dict[c] for c in self._wb_header]))
        self._t = self.get_current_time()

    def update_until(self, t, dt=None):
        """Updates the model until time t is reached with an update time step dt.
        
        Arguments:
            t {int} -- Model time until which the model is update
        
        Keyword Arguments:
            dt {int} -- update time step interval [s] at which information is exchanged between models (default: {None})
        
        Raises:
            Exception -- Raised if specified time is smaller than time step or later than model end time
        """

        if (t<self._t) or t>self.get_end_time():
            msg = "wrong time input: smaller than model time or larger than endTime"
            self.logger.error(msg)
            raise Exception(msg)
        while self._t < t:
            self.update(dt=dt)

    def spinup(self):
        """PCR-GLOBWB specific function (see pcr_bmi.py).
        
        Raises:
            NotImplementedError -- SpinUp is not implemented for central BMI
        """

        raise NotImplementedError

    def finalize(self):
        """Finalizes the model.
        shuts down all operations and closes output files.
        """

        self.logger.info('finalize models bmi. Close loggers.')
        for mod in self.bmimodels:
            self.bmimodels[mod].finalize()
        closelogger(self.logger)
        closelogger(self.wb_logger)

    def exchange(self, from_mod, to_mod, from_vars, to_vars, coupling, add=False, **kwargs):
        """Exchanges variable content from specified source variable in source model to
        a specified destination variable in the destination model. 
        exchange is possible for entire grid or at specified indices;
        coupling either as totals or fractional values depending on exchange;
        source variable can either add to destination variable or overwrite it
        
        Arguments:
            from_mod {str} -- string defining the source model
            to_mod {str} -- string defining the destination model
            from_vars {str} -- string defining the source variable
            to_vars {str} -- string defining the destination variable
            coupling {SpatialCoupling} -- object defining the spatial coupling structure
        
        Keyword Arguments:
            add {bool} -- if True, source values are added to destination values; if False, overwritten (default: {False})
        """

        self.logger.info('{} {}.{} data to coupled {}.{} variable'.format('adding' if add else 'setting', from_mod, from_vars[0], to_mod, to_vars[0]))
        if coupling.at_indices():
            tot_volume = self.exchange_at_indices(from_mod, to_mod, from_vars, to_vars, coupling, add=add, **kwargs)
        else:
            tot_volume = self.exchange_same_grid(from_mod, to_mod, from_vars, to_vars, coupling, add=add, **kwargs)
        return tot_volume

    def exchange_same_grid(self, from_mod, to_mod, from_vars, to_vars, coupling, add=False, **kwargs):
        """Exchanges variable content from specified source variable in source model to
        a specified destination variable in the destination model for the entire grid.
        source variable can either add to destination variable or overwrite it
        
        Arguments:
            from_mod {str} -- string defining the source model
            to_mod {str} -- string defining the destination model
            from_vars {str} -- string defining the source variable
            to_vars {str} -- string defining the destination variable
            coupling {SpatialCoupling} -- object defining the spatial coupling structure
        Keyword Arguments:
            add {bool} -- if True, source values are added to destination values; if False, overwritten (default: {False})
        """

        # get FROM data & translate
        shape = self.bmimodels[from_mod].get_var_shape(from_vars[0])
        vals = np.ones(shape)
        for from_var in from_vars:
            if isinstance(from_var, float): # scalar multiplier
                vals *= from_var
            else:
                vals *= self.bmimodels[from_mod].get_value(from_var)
        # reproject vals to the to_mod projection using the SpatialCoupling.reproject function
        if coupling.reproject is not None:
            vals = coupling.reproject(vals, np.nan)
        # total exchange
        self.logger.debug('total get {:3f}'.format(np.nansum(vals)))
        tot_volume = np.nansum(vals)
        # this is an ugly hack at the moment to circumvent the PCR.runoff [m] to CMF.roffin [m] coupling
        # TODO: solve this when we implement pint for unit conversion
        if (from_mod == 'PCR') and (from_vars[0] == 'runoff') and (to_mod == 'CMF') and (to_vars[0] == 'roffin'):
            tot_volume = np.nansum(vals * self.bmimodels[from_mod].get_value(self.bmimodels[from_mod]._area_var_name))
        # get TO data & translate
        for var in to_vars[1:]:
            if isinstance(var, float):
                div = var
            else:
                div = self.bmimodels[to_mod].get_value(var) 
            assert np.all(np.not_equal(div, 0)) # make sure not to divide by zero
            vals /= div
        # SET data
        # ensure no missing values appear
        vals[np.isnan(vals)] = 0.
        vals = np.maximum(vals, 0)
        if add:  # add to current state
            vals += self.bmimodels[to_mod].get_value(to_vars[0])

        self.bmimodels[to_mod].set_value(to_vars[0], vals)
        # import matplotlib.pyplot as plt
        # import pdb;pdb.set_trace()
        return tot_volume

    def exchange_at_indices(self, from_mod, to_mod, from_vars, to_vars, coupling, add=False, **kwargs):
        """Exchanges variable content from specified source variable in source model to
        a specified destination variable in the destination model for user-specified indices.
        coupling either as totals or fractional values depending on exchange;
        source variable can either add to destination variable or overwrite it
           
        Arguments:
            from_mod {str} -- string defining the source model
            to_mod {str} -- string defining the destination model
            from_vars {str} -- string defining the source variable
            to_vars {str} -- string defining the destination variable
            coupling {SpatialCoupling} -- object defining the spatial coupling structure
        
        Keyword Arguments:
            add {bool} -- if True, source values are added to destination values; if False, overwritten (default: {False})
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
        tot_volume = np.nansum(vals_get)
        vol_diff = tot_volume - np.nansum(vals_set)
        if vol_diff > 1E-2:
            self.logger.warming('Large difference in water volume from and to: {:.4f} m3'.format(vol_diff))
        self.logger.debug('total get {:3f}'.format(tot_volume))
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
        return tot_volume

    ###
    ### Model Information Functions
    ###

    def get_model_type(self):
        """Returns type of model. 
        possible types are hydrologic, routing, or hydrodynamic model
        
        Returns:
            str -- type of specified model
        """

        return {mod: self.bmimodels[mod].get_model_type() for mod in self.bmimodels}

    def get_component_name(self):
        #TODO: what are model component in this context?
        """Returns component name of specified model.
        
        Returns:
            str -- name of specified component
        """

        return {mod: self.bmimodels[mod].get_component_name() for mod in self.bmimodels}

    def get_input_var_names(self):
        """Returns list with all possible input variable names that 
        can be used as exchange TO the specified model.
        must be specified in model specific BMI class.
        
        Returns:
            list -- list of all possible input variables
        """

        return {mod: self.bmimodels[mod].get_input_var_names() for mod in self.bmimodels}

    def get_output_var_names(self):
        """Returns list with all possible output variable names that 
        can be used as exchange FROM the specified model.
        must be specified in model specific BMI class.
        
        Returns:
            list -- list of all possible output variables
        """

        return {mod: self.bmimodels[mod].get_output_var_names() for mod in self.bmimodels}

    ###
    ### Variable Information Functions
    ###

    def get_var_type(self, long_var_name):
        """Returns the type of a user-specified model variable exposed via BMI.
        
        Arguments:
            long_var_name {str} -- long name of exposed model variable
        
        Returns:
            str -- type of specified variable
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_type(var)

    def get_var_units(self, long_var_name):
        """Provides units of variable.
        
        Arguments:
            long_var_name {str} -- long name of exposed model variable
        
        Returns:
            str -- units of specified variable
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_units(var)

    def get_var_rank(self, long_var_name):
        """Provides number of dimensions of variable.
        
        Arguments:
            long_var_name {str} -- long name of exposed model variable
        
        Returns:
            int -- rank of specified variable
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_rank(var)

    def get_var_size(self, long_var_name):
        """Providestotal number of values contained in variable.
        
        Arguments:
            long_var_name {str} -- long name of exposed model variable
        
        Returns:
            int -- size of specified variable
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_size(var)

    def get_var_shape(self, long_var_name):
        """Provides shape of variable.
        
        Arguments:
            long_var_name {str} -- long name of exposed model variable
        
        Returns:
            tuple -- shape of specified variable
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_shape(var)

    def get_var_nbytes(self, long_var_name):
        """Provides number of bytes of variable.
        
        Arguments:
            long_var_name {str} -- long name of exposed model variable
        
        Returns:
            int -- byte number of specified variable
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_var_nbytes(var)

    def get_start_time(self):
        """Provides start time of model.
        
        Returns:
            date -- start time of model
        """

        self._startTime = max([self.bmimodels[mod].get_start_time() for mod in self.bmimodels])
        return self._startTime

    def get_current_time(self):
        """Provides current model time.
        
        Returns:
            date -- current model time of PCR-GLOBWB
        """

        models_t = [self.bmimodels[mod]._t for mod in self.bmimodels]
        in_sync = not models_t or models_t.count(models_t[0]) == len(models_t)  # check if identical
        if not in_sync:
            msg = '; '.join(['{}: {}'.format(m, str(t)) for m, t in zip(self.bmimodels, models_t)])
            self.logger.warning("Model times out of sync: {}".format(msg))
        self._t = models_t[0]
        return self._t

    def get_end_time(self):
        """Provides end time of model.
        
        Returns:
            date -- end time of model
        """

        self._endTime = min([self.bmimodels[mod].get_end_time() for mod in self.bmimodels])
        return self._endTime

    def get_time_step(self):
        """Provides time step of model.
        
        Returns:
            date -- time step of model
        """
        
        return self._dt
        
    def get_time_units(self):
        """Provides time unit of model.
        
        Returns:
            str -- time unit of model
        """

        return self._timeunit

    ###
    ### Variable Getter and Setter Functions
    ###
    
    def get_value(self, long_var_name, **kwargs):
        """Gets values of a certain exposed variable. 
        entries with fill_value are replaced with NaN values
        
        Arguments:
            long_var_name {str} -- name of exposed model variable
        
        Returns:
            array -- array with values
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_value(var, **kwargs)

    def get_value_at_indices(self, long_var_name, inds, **kwargs):
        """Gets values at a specific index of a certain exposed variable. 
        entries with fill_value are replaced with NaN values
        
        Arguments:
            long_var_name {str} -- name of exposed model variable
            inds {int} -- Index pointing to entry within entire array of variable values
        
        Returns:
            float -- value at specific index of retrieved variable
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].get_value_at_indices(var, inds, **kwargs)

    def set_value(self, long_var_name, src, **kwargs):
        """Overwriting of all values of a certain exposed variable with provided new values. 
        entries with NaN value are replaced with fill_value; 
        provided new values must match shape of aim variable
        
        Arguments:
            long_var_name {str} -- name of exposed model variable
            src {array} -- array with new values
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].set_value(var, src, **kwargs)

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        """Overwriting of value at specific entry of a certain exposed variable with provided new values. 
        entries with NaN value are replaced with fill_value
        
        Arguments:
            long_var_name {str} -- name of exposed model variable
            inds {int} -- Index pointing to entry within entire array of variable values
            src {array} -- array with new values
        
        Returns:
            [type] -- [description]
        """

        mod, var = self._check_long_var_name(long_var_name)
        return self.bmimodels[mod].set_value_at_indices(var, inds, src, **kwargs)

    ###
    ### Grid Information Functions
    ###
    
    def get_grid_type(self):
        """Provides grid type of model.
        can be regular, flexible or unit catchment grid
        
        Returns:
            str -- Grid type of model
        """

        return {mod: self.bmimodels[mod].get_grid_type() for mod in self.bmimodels}

    ###
    ### Observation points
    ###

    def get_obs(self):
        """Under construction;
        function to get observation points to read out model
        output on the fly.
        """

        pass

    def set_obs(self):
        """Under construction;
        function to set observation points to read out model
        output on the fly.
        """

        pass

    def index(self, x, y, mod, in1d=False):
        """Finds the index in model data grid associated with
        coordinates in lat/lon projection of a user-specified point.
        
        Arguments:
            x {float} -- x-coordinate of point
            y {float} -- y-coordinate of point
            mod {str} -- model abbreviation
        
        Keyword Arguments:
            in1d {bool} -- True if point is on a 1D network in model (default: {False})
        
        Raises:
            AssertionError -- Raised if no grid is set before
            AssertionError -- Raised if no 1D network is set before
        
        Returns:
            int -- index associated to x/y coordinates
        """

        if self.bmimodels[mod].grid is None:
            msg = "{}: model grid not set".format(mod)
            self.logger.error(msg)
            raise AssertionError(msg)
        if in1d:
            if self.bmimodels[mod].grid._1d is None:
                msg = "{}: model 1d network not set".format(mod)
                self.logger.error(msg)
                raise AssertionError(msg)
            idx = self.bmimodels[mod].grid._1d.index(x,y)
        else:
            idx = self.bmimodels[mod].grid.index(x,y)
        return idx
    
    ###
    ### set and get attribute / config 
    ###

    def set_start_time(self, start_time):
        """Overwriting default model start time with user-specified start time.
        format of provided start time must be yyyy-mm-dd
        
        Arguments:
            start_time {date} -- user-specified start time
        """

        for mod in self.bmimodels:
            self.bmimodels[mod].set_start_time(start_time)
        self.get_start_time() # sync start time
        self._t = self._startTime
       
    def set_end_time(self, end_time):
        """Overwriting default model end time with user-specified end time.
        format of provided end time must be yyyy-mm-dd
        
        Arguments:
            end_time {date} -- user-specified end time
        """

        for mod in self.bmimodels:
            self.bmimodels[mod].set_end_time(end_time)
        self.get_end_time() # sync end time

    def get_attribute_names(self):
        """Provides list with all model attribute names from model config file.
        
        Returns:
            list -- list with model attribute names
        """

        return {mod: self.bmimodels[mod].get_attribute_names() for mod in self.bmimodels}

    def set_out_dir(self, out_dir):
        """Setting output directory of model. 
        overwrites the default output directory
        
        Arguments:
            out_dir {str} -- path to output directory
        """

        for mod in self.bmimodels:
            path = join(out_dir, mod)
            if not os.path.isdir(path):
                os.makedirs(path)
            self.bmimodels[mod].set_out_dir(path)

    def get_attribute_value(self, attribute_name):
        """gets attribute value in in underlying model.
        attribute_name should use the following convention model_name.section_name:attribute_name
        
        Arguments:
            attribute_name {str} -- Name of BMI attribute
        
        Returns:
            float -- value of BMI attribute
        """

        mod, attr = self._check_long_var_name(attribute_name)
        return self.bmimodels[mod].get_attribute_value(attr)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        """sets attribute value in underlying model.
        attribute_name should use the following convention model_name.section_name:attribute_name
        
        Arguments:
            attribute_name {str} -- Name of BMI attribute
            attribute_value {float} -- value to be set
        """

        mod, attr = self._check_long_var_name(attribute_name)
        return self.bmimodels[mod].set_attribute_value(attr, attribute_value)

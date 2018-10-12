import numpy as np
import sys
from configparser import ConfigParser
import logging
import os
from os.path import join, isfile, abspath, dirname, basename, normpath
from datetime import datetime, timedelta
import rasterio

from utils import setlogger, closelogger
from gbmi import GBmi
from grids import RGrid
import glofrim_lib as glib 


class PCR(GBmi):
    """CSDMS-compliant BMI implementation for GLOFRIM.
    
    Arguments:
        GBmi {class} -- Interface (abstract base class) for a model that implements the CSDMS BMI (Basic Model Interface).
    
    Raises:
        Warning -- Warning is raised if two-step model initialization is not correctly executed
        Warning -- Warning is raised if two-step model initialization is not correctly executed
        Exception -- Raised if specified time is smaller than time step or later than model end time
        ValueError -- Raised if model time step is not whole number of update interval
        Exception -- Raised if end time is already reached as no further update is possible
        IOError -- function requires LDD file; if not found, IOerror is raised
        IOError -- function requires landsmask file; if not found, IOerror is raised
        ValueError -- start time must be yyyy-mm-dd; if not, ValueError is raised
        ValueError -- end time must be yyyy-mm-dd; if not, ValueError is raised
    """

    _name = 'PCR'
    _long_name = 'PCR-GLOBWB'
    _version = '2.0.3'
    _var_units = {'runoff': 'm/day', 'discharge': 'm3/s', 'cellArea': 'm2'}
    _input_var_names = ['cellArea']
    _output_var_names = ['runoff', 'discharge']
    _area_var_name = 'cellArea'
    _timeunit = 'days'
    _dt = timedelta(days=1) # NOTE: this is fixed in PCR

    def __init__(self, loglevel=logging.INFO, logger=None):
        # import PCR-GLOBWB with BMI functions 
        from pcrglobwb_bmi_v203 import pcrglobwb_bmi as _bmi
        self._bmi = _bmi.pcrglobwbBMI()
        if logger:
            self.logger = logger.getChild(self._name)
        else:
            self.logger = setlogger(None, self._name, thelevel=loglevel)
        self._loglevel = loglevel
        self.initialized = False
        self.grid = None

    ###
    ### Model Control Functions
    ###

    def initialize_config(self, config_fn):
        """Initializing the model configuration file. Aligning GLOFRIM specs for coupled runs
        to be consistent with overall run settings.
        
        Arguments:
            config_fn {str} -- path to model configuration file (for PCR-GLOBWB: ini-file)
        
        Raises:
            Warning -- Warning is raised if two-step model initialization is not correctly executed
        """

        # config settings
        if self.initialized:
            raise Warning("model already initialized, it's therefore not longer possible to initialize the config")
        self._config_fn = abspath(config_fn)
        self._config = glib.configread(self._config_fn, encoding='utf-8', 
                                       cf=ConfigParser(inline_comment_prefixes=('#')))
        self._datefmt = "%Y-%m-%d"
        # model time
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model files
        self._outdir= abspath(self.get_attribute_value('globalOptions:outputDir'))
        self.logger.info('Config initialized')

    def initialize_model(self, **kwargs):
        """Initializes the model, i.e. loading all files and checking for consistency; 
        can only be executed after the model-specific config-file was initialized. This is essential
        for a successful two-step model initialization and aligned model data.
        
        Raises:
            Warning -- Warning is raised if two-step model initialization is not correctly executed
        """

        if not hasattr(self, '_config_fn'):
            raise Warning('run initialize_config before initialize_model')
        self.write_config() # write updated config to file as bmi does not allow direct access
        self._bmi.initialize(self._config_fn)
        # stop pcr double logging.
        pcr_log = logging.getLogger()
        # pcr_log.setLevel(self._loglevel)
        self.logger.info("remove PCR logger because it's making to much noise")
        pcr_log.handlers = pcr_log.handlers[1:]
        self.initialized = True
        self.logger.info('Model initialized')
        # reset model time to make sure it is consistent with the model
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime

    def initialize(self, config_fn):
        """Initializes the model following a two-step initialization procedure.
        First, the config-file is initialized and where necessary aligned with
        overall model settings (e.g. output-dir, start and end time).
        Second, the model is actually initialized.
        
        Arguments:
            config_fn {str} -- path to model configuration file (for PCR-GLOBWB: ini-file)
        """

        self.initialize_config(config_fn)
        self.initialize_model()

    def update(self, dt=None):
        """Updating model for a certain time step interval (default: None).
        Checks whether model end time is already reached.
        Requires that model-specific time step is a whole number of update time step.
        
        Keyword Arguments:
            dt {int} -- update time step interval [s] at which information is exchanged between models (default: {None})
        
        Raises:
            ValueError -- Raised if model time step is not whole number of update interval
            Exception -- Raised if end time is already reached as no further update is possible
        """

        # dt in seconds. if not given model timestep is used
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        if (dt is not None) and (dt != self._dt.total_seconds()):
            _dt = timedelta(seconds=dt)
            if not glib.check_dts_divmod(_dt, self._dt):
                msg = "Invalid value for dt in comparison to model dt. Make sure a whole number of model timesteps ({}) fit in the given dt ({})"
                raise ValueError(msg.format(self._dt, _dt))
            self._bmi.update(_dt.days)
            self._t += _dt
        else:
            self._bmi.update(self._dt.days)
            self._t += self._dt
        self.logger.info('updated model to datetime {}'.format(self._t.strftime("%Y-%m-%d %H:%M")))

    def update_until(self, t, dt=None):
        """Updates the model until a certain time t is reached with an update time step dt.
        
        Arguments:
            t {int} -- Model time until which the model is updated
        
        Keyword Arguments:
            dt {int} -- update time step interval [s] at which information is exchanged between models (default: {None})
        
        Raises:
            Exception -- Raised if specified time is smaller than time step or later than model end time
        """

        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update(dt=dt)

    def spinup(self):
        """PCR-GLOBWB specific function.
        Runs the in ini-file specified number of spin-up years.
        """

        self._bmi.spinup()

    def finalize(self):
        """Finalizes the model, i.e. shuts down all operations and closes output files.
        """

        self.logger.info('finalize bmi. Close logger.')
        self._bmi.finalize()
        closelogger(self.logger)

    ###
    ### Variable Information Functions
    ###
    
    def get_start_time(self):
        """Provides start time of PCR-GLOBWB
        
        Returns:
            date -- start time of PCR-GLOBWB
        """

        if self.initialized:
            # date to datetime object
            startTime = datetime.combine(self._bmi.get_start_time(), datetime.min.time())
        else:
            startTime = self.get_attribute_value('globalOptions:startTime')
            startTime = datetime.strptime(startTime, self._datefmt)
        self._startTime = startTime
        return self._startTime
    
    def get_current_time(self):
        """Provides current model time of PCR-GLOBWB
        
        Returns:
            [date] -- current model time of PCR-GLOBWB
        """

        if self.initialized:
            return self._t
        else:
            return self.get_start_time()

    def get_end_time(self):
        """Provides end time of PCR-GLOBWB
        
        Returns:
            date -- end time of PCR-GLOBWB
        """

        if self.initialized:
            # date to datetime object
            endTime = datetime.combine(self._bmi.get_end_time(), datetime.min.time())
        else:
            endTime = self.get_attribute_value('globalOptions:endTime')
            endTime = datetime.strptime(endTime, self._datefmt)
        self._endTime = endTime
        return self._endTime

    def get_time_step(self):
        """Provides time step of PCR-GLOBWB
        
        Returns:
            date -- time step of PCR-GLOBWB
        """

        return self._dt 
        
    def get_time_units(self):
        """Provides time unit of PCR-GLOBWB.
        
        Returns:
            str -- time unit of model of PCR-GLOBWB
        """

        return self._timeunit

    ###
    ### Variable Getter and Setter Functions
    ###
    
    def get_value(self, long_var_name, fill_value=-999, **kwargs):
        """Retrieval of all values of a certain exposed variable; 
        entries with fill_value are replaced with NaN values
        
        Arguments:
            long_var_name {str} -- name of exposed PCR-GLOBWB variable
        
        Keyword Arguments:
            fill_value {int} -- fill_value used in PCR-GLOBWB (default: {-999})
        
        Returns:
            array -- array with all entries of retrieved variable
        """

        # additional fill_value argument required to translate pcr maps to numpy arrays
        var = np.asarray(self._bmi.get_var(long_var_name, missingValues=fill_value, **kwargs))
        var = np.where(var == fill_value, np.nan, var)
        return var

    def get_value_at_indices(self, long_var_name, inds, fill_value=-999, **kwargs):
        """Retrieval of value at a specific index of a certain exposed variable; 
        entries with fill_value are replaced with NaN values
        
        Arguments:
            long_var_name {str} -- name of exposed PCR-GLOBWB variable
            inds {int} -- Index pointing to entry within entire array of variable values
        
        Keyword Arguments:
            fill_value {int} -- fill_value used in PCR-GLOBWB (default: {-999})
        
        Returns:
            float -- value at specific index of retrieved variable
        """

        return self.get_value(long_var_name, fill_value=fill_value, **kwargs).flat[inds]

    def set_value(self, long_var_name, src, fill_value=-999, **kwargs):
        """Overwriting of all values of a certain exposed variable with provided new values; 
        entries with NaN value are replaced with fill_value; 
        provided new values must match shape of aim variable
        
        Arguments:
            long_var_name {str} -- name of exposed PCR-GLOBWB variable
            src {array} -- array with new values
        
        Keyword Arguments:
            fill_value {int} -- fill_value used in PCR-GLOBWB (default: {-999})
        """

        src = np.where(np.isnan(src), fill_value, src).astype(self.get_var_type(long_var_name))
        self._bmi.set_var(long_var_name, src, missingValues=fill_value, **kwargs)

    def set_value_at_indices(self, long_var_name, inds, src, fill_value=-999,  **kwargs):
        """Overwriting of value at specific entry of a certain exposed variable with provided new values; 
        entries with NaN value are replaced with fill_value
        
        Arguments:
            long_var_name {str} -- name of exposed PCR-GLOBWB variable
            inds {int} -- Index pointing to entry within entire array of variable values
            src {array} -- array with new values
        
        Keyword Arguments:
            fill_value {int} -- fill_value used in PCR-GLOBWB (default: {-999})
        """

        val = self.get_value(long_var_name, fill_value=fill_value, **kwargs)
        val.flat[inds] = src
        self.set_value(long_var_name, val, **kwargs)
        
    ###
    ### Grid Information Functions
    ###

    def get_grid(self):
        """Retreving spatial information about model domain from PCR-GLOBWB landmask file using rasterio
        operations; 
        also retrieving PCR-GLOBWB drainage direction (LDD) information
        
        Raises:
            IOError -- function requires LDD file; if not found, IOerror is raised
            IOError -- function requires landsmask file; if not found, IOerror is raised
        
        Returns:
            grid -- rasterio grid object
        """


        if not hasattr(self, 'grid') or (self.grid is None):
            # ldd is used for coupling to routing / flood model
            _indir = abspath(self.get_attribute_value('globalOptions:inputDir'))
            _ldd_fn = glib.getabspath(self.get_attribute_value('routingOptions:lddMap'), _indir)
            if not isfile(_ldd_fn): raise IOError('ldd file not found {}'.format(_ldd_fn))
            # landmask used for masking out coordinates outside mask
            _lm_fn = glib.getabspath(self.get_attribute_value('globalOptions:landmask'), _indir)
            if isfile(_lm_fn): #raise IOError('landmask file not found')
                with rasterio.open(_lm_fn, 'r') as ds:
                    mask=ds.read(1)==1
                    mask_name = basename(_lm_fn)
            else:
                mask=None
                mask_name = basename(_ldd_fn)
            self.logger.info('Getting rgrid info based on {}; mask based on {}'.format(basename(_ldd_fn), mask_name))
            with rasterio.open(_ldd_fn, 'r') as ds:
                if mask is None:
                    mask=ds.read(1) != ds.nodata
                self.grid = RGrid(ds.transform, ds.height, ds.width, crs=ds.crs, mask=mask)
            # read file with pcr readmap
            self.logger.info('Getting drainage direction from {}'.format(basename(_ldd_fn)))
            self.grid.set_dd(_ldd_fn, ddtype='ldd')
        return self.grid

    ###
    ### set and get attribute / config 
    ###

    def set_start_time(self, start_time):
        """Overwriting default model start time with user-specified start time;
        format of provided start time must be yyyy-mm-dd
        
        Arguments:
            start_time {date} -- user-specified start time
        
        Raises:
            ValueError -- start time must be yyyy-mm-dd; if not, ValueError is raised
        """

        if isinstance(start_time, datetime):
            start_time = start_time.strftime(self._datefmt)
        try:
            self._startTime = datetime.strptime(start_time, self._datefmt) 
            self._t = self._startTime
        except ValueError:
            raise ValueError('wrong date format, use "yyyy-mm-dd"')
        self.set_attribute_value('globalOptions:startTime', start_time)
       
    def set_end_time(self, end_time):
        """Overwriting default model end time with user-specified end time;
        format of provided end time must be yyyy-mm-dd
        
        Arguments:
            end_time {date} -- user-specified end time
        
        Raises:
            ValueError -- end time must be yyyy-mm-dd; if not, ValueError is raised
        """

        if isinstance(end_time, datetime):
            end_time = end_time.strftime(self._datefmt)
        try:
            self._endTime = datetime.strptime(end_time, self._datefmt) 
        except ValueError:
            raise ValueError('wrong date format, use "yyyy-mm-dd"')
        self.set_attribute_value('globalOptions:endTime', end_time)

    def set_out_dir(self, out_dir):
        """Setting output directory of PCR-GLOBWB; 
        overwrites the default output directory
        
        Arguments:
            out_dir {str} -- path to output directory
        """

        self.set_attribute_value('globalOptions:outputDir', abspath(out_dir))

    def get_attribute_names(self):
        """Provides list with all model attribute names from model config file
        
        Returns:
            list -- list with model attribute names
        """

        glib.configcheck(self, self.logger)
        return glib.configattr(self._config)
    
    def get_attribute_value(self, attribute_name):
        """Gets attribute value in in underlying model;
        attribute_name should use the following convention model_name.section_name:attribute_name
        
        Arguments:
            attribute_name {str} -- Name of BMI attribute
        
        Returns:
            float -- value of attribute
        """

        glib.configcheck(self, self.logger)
        self.logger.debug("get_attribute_value: " + attribute_name)
        return glib.configget(self._config, attribute_name)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        """sets attribute value in in underlying model;
        attribute_name should use the following convention model_name.section_name:attribute_name
        
        Arguments:
            attribute_name {str} -- Name of BMI attribute
            attribute_value {float} -- value to be est
        
        Returns:
            [str] -- [no clue]
        """

        glib.configcheck(self, self.logger)
        self.logger.debug("set_attribute_value: " + attribute_value)
        return glib.configset(self._config, attribute_name, attribute_value)

    def write_config(self):
        """Write adapted config to file; just before initializing;
        only for models which do not allow for direct access to model config via bmi
        """

        self._config_fn = glib.write_config(self, self._config, self._config_fn, self.logger)
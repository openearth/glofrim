# -*- coding: utf-8 -*-
import numpy as np
from pcrglobwb_bmi_v203 import pcrglobwb_bmi
from bmi.wrapper import BMIWrapper
import logging
import warnings

logger = logging.getLogger(__name__)

# wrapper around BMI
class _model(object):
    def __init__(self, bmi, name, mtype, missing_value=np.nan):
        self.bmi = bmi
        self.name = name
        self.mtype = mtype
        self.missing_value = missing_value
        self.start_time = self.get_start_time()

    def get_var(self, var, parse_missings=True, *args, **kwargs):
        var_data = self.bmi.get_var(var)
        if parse_missings:
            # if given nodata is parsed to np.nan
            var_data = np.where(var_data == self.missing_value, np.nan, var_data)
        return var_data

    def set_var(self, *args, **kwargs):
        self.bmi.set_var(*args, **kwargs)

    def spinup(self, *args, **kwargs):
        self.bmi.spinup(*args, **kwargs)

    def get_start_time():
        "get model start (initialization) time"
        self.bmi.get_start_time()

    def get_current_time():
        "get model current time"
        self.bmi.get_start_time()

    def get_time_step():
        "get model current time step"
        self.bmi.get_start_time()

    def update(self, dt=-1):
        self.bmi.update(dt=dt)
        start_time = self.start_time
        current_time = self.get_current_time()
        time_step = self.bmi.get_time_step()
        logger.info(
            "%s -> start_time: %s, current_time %s, timestep %s",
            self.name,
            start_time, #TODO replace by self.start_time
            current_time,
            time_step
        )


class PCR_model(_model):
    def __init__(self, configfile, missing_value):
        pcr_bmi = pcrglobwb_bmi.pcrglobwbBMI()
        pcr_bmi.initialize(configfile)
        super(PCR_model, self).__init__(pcr_bmi, 'PCRGLOB-WB', 'hydrology',
                                        missing_value=missing_value)

    def run(self, dt=1):
        "run model for dt [days]"
        super(PCR_model, self).update(dt=dt)

class CMF_model(_model):
    def __init__(self, engine, configfile, map_dir, get_runoff=None):
        """GLOFRIM wrapper around BMI for CaMa-Flood (CMF)

        input
        engine          BMIWrapper path to CMF engine
        configfile      BMIWrapper CMF config file name
        map_dir         CMF map directory. Required to initialize CMF
        get_runoff      Function to get runoff from upstream coupled model.
                        This can also be set later using the set_options function.

        """
        cmf_bmi = BMIWrapper(engine = engine, configfile = configfile)
        cmf_bmi.initialize(map_dir)
        super(CMF_model, self).__init__(cmf_bmi, 'CaMa-Flood', 'hydrodynamic')
        self.get_runoff = get_runoff

    def set_options(self, get_runoff=None):
        if get_runoff is not None:
            self.get_runoff = get_runoff
            logger.info('CMF get runoff function set')


    def run(self, dt=1):
        "update runoff state and run for dt [days]"
        dt = dt * 86400 # convert days to seconds
        if self.get_runoff is not None:
            self.set_var('runoff',  self.get_runoff())
            super(CMF_model, self).update(dt=dt)
        else:
            warnings.warn(
                "set get_runoff function first using the set_options function",
                RuntimeWarning
            )

    def get_var(self, var, parse_missings=True, *args, **kwargs):
        var = super(CMF_model, self).get_var(var, parse_missings=parse_missings,
                                             *args, **kwargs)
        # return var with switched axis (fortran to python translation)
        return var.reshape(var.shape[::-1])

    def _CMFtime_2_datetime(self, t, yr0=1850):
        # internal CMF time is in minutes since 1850
        h, m = divmod(t, 60)
        d, h = divmod(h, 24)
        return datetime(yr0, 1, 1) + timedelta(days = d)

    def get_current_time(self):
        t = super(CMF_model, self).get_current_time()
        return self._CMFtime_2_datetime(t)

    def get_start_time(self):
        t = super(CMF_model, self).get_start_time()
        return self._CMFtime_2_datetime(t)

# class DFM_model(_model):
#     def __init__(self, path_to_model, configfile):
#         cmf_bmi = BMIWrapper(engine = path_to_model, configfile = configfile)
#         cmf_bmi.initialize()
#         super(CMF_model, self).__init__(cmf_bmi, 'Delft3D-FM', 'hydrodynamic')
#
#     def update(self, runoff, *args, **kwargs):
#         self.set_var('runoff', runoff)
#         self.bmi.update(*args, **kwargs)

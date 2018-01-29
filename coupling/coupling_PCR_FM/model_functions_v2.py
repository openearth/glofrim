# -*- coding: utf-8 -*-
import numpy as np
from pcrglobwb_bmi_v203 import pcrglobwb_bmi
from bmi.wrapper import BMIWrapper
import logging
import warnings
from datetime import datetime, timedelta
import os, shutil
from utils import write_ini
logger = logging.getLogger(__name__)

# wrapper around BMI
class _model(object):
    def __init__(self, bmi, name, missing_value=np.nan):
        self.bmi = bmi
        self.name = name
        self.missing_value = missing_value
        self.start_time = self.get_start_time()
        #

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

    def get_start_time(self):
        "get model start (initialization) time"
        return self.bmi.get_start_time()

    def get_current_time(self):
        "get model current time"
        return self.bmi.get_current_time()

    def get_time_step(self):
        "get model current time step"
        return self.bmi.get_time_step()

    def run(self, dt=1):
        "run model for dt [days]"
        self.update(dt=dt)

    def update(self, dt=-1):
        self.bmi.update(dt=dt)
        current_time = self.get_current_time()
        time_step = self.get_time_step()
        logger.info(
            "%s -> start_time: %s, current_time %s, timestep %s",
            self.name,
            self.start_time, #TODO replace by self.start_time
            current_time,
            time_step
        )


class PCR_model(_model):
    def __init__(self, config_fn,
                 start_date, end_date,
                 out_dir, in_dir,
                 missing_value=-999,
                 ):
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        # make ini file consistent with GLOFRIM settings
        date_fmt = "%Y-%m-%d"
        ini_kwargs = {"inputDir": os.path.abspath(in_dir),
                      "outputDir": os.path.abspath(out_dir),
                      "startTime": start_date.strftime(date_fmt),
                      "endTime": end_date.strftime(date_fmt),
                      }
        tmp_config_fn = os.path.join(out_dir, 'tmp_' + os.path.basename(config_fn))
        write_ini(tmp_config_fn, config_fn, **ini_kwargs)
        logger.info('tmp ini file for PCR model written to {:s}'.format(tmp_config_fn))

        # initialize BMIWrapper and model
        pcr_bmi = pcrglobwb_bmi.pcrglobwbBMI()
        # TODO: check why pcr_bmi changes cwd. it couses errors downstream
        pcr_bmi.initialize(tmp_config_fn)
        super(PCR_model, self).__init__(pcr_bmi, 'PCRGLOB-WB',
                                        missing_value=missing_value)


class CMF_model(_model):
    def __init__(self, engine, config_fn, model_dir,
                 start_date, end_date,
                 out_dir, in_dir=None,
                 missing_value=1e20):
        """GLOFRIM wrapper around BMI for CaMa-Flood (CMF)

        input
        engine          BMIWrapper path to CMF engine
        config_fn      BMIWrapper CMF config file name
        map_dir         CMF map directory. Required to initialize CMF
        update_states   Function to get update the states (e.g. runoff) based on coupled model.
                        This can also be set later using the set_options function.

        """
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        ## make inifile consisitent with settings
        ini_kwargs = {"COUTDIR": os.path.abspath(out_dir), # output dir
                      "ISYEAR": "{:d}".format(start_date.year),  # model start & end dates
                      "ISMON": "{:d}".format(start_date.month),
                      "ISDAY": "{:d}".format(start_date.day),
                      "IEYEAR": "{:d}".format(end_date.year),
                      "IEMON": "{:d}".format(end_date.month),
                      "IEDAY": "{:d}".format(end_date.day),
                      }
        if in_dir is None: # do not read from file but BMI
            # TODO: write consistent inpmat file and add to write_ini
            ini_kwargs.update({"CRUNOFFDIR": "",
                               "LBMIROF": ".TRUE.",
                               "DTIN": "86400", # input unit settings -> runoff in m/s
                               "DROFUNIT": "1",
                               # "CDIMINFO":  # input dimensions file
                               # "CINPMAT":
                               })
        tmp_config_fn = os.path.join(model_dir, "input_flood.nam")
        if os.path.abspath(tmp_config_fn) == os.path.abspath(config_fn):
            shutil.move(config_fn, config_fn + '.template')
            config_fn = config_fn + '.template'
            logger.warning('template .nam file cannot be called input_flood.nam; file renamed {:s}'.format(config_fn))
        write_ini(tmp_config_fn, config_fn, **ini_kwargs)
        logger.info('tmp ini file for CMF model written to {:s}'.format(tmp_config_fn))

        ## initialize BMIWrapper and model
        # NOTE that that the CMF BMI wrapper does not use the
        cmf_bmi = BMIWrapper(engine = engine)
        cmf_bmi.initialize(model_dir)
        super(CMF_model, self).__init__(cmf_bmi, 'CaMa-Flood',
                                        missing_value=missing_value)
        self.update_states = None

    def set_options(self, update_states=None):
        if update_states is not None:
            self.update_states = update_states
            logger.info('CMF update_states function set')


    def run(self, dt=1):
        "update CMF states and run for dt [days]"
        dt = dt * 86400 # convert days to seconds
        if self.update_states is not None:
            self.update_states()
        else:
            warnings.warn(
                "Update_states function not set. No states updated before runtime",
                RuntimeWarning
            )
        super(CMF_model, self).update(dt=dt)

    def get_var(self, var, parse_missings=True, *args, **kwargs):
        var = super(CMF_model, self).get_var(var, parse_missings=parse_missings,
                                             *args, **kwargs)
        # return var with switched axis (fortran to python translation)
        return var.reshape(var.shape[::-1])

    def get_current_time(self):
        t = super(CMF_model, self).get_current_time()
        return CMFtime_2_datetime(t)

    def get_start_time(self):
        t = super(CMF_model, self).get_start_time()
        return CMFtime_2_datetime(t)


# utils
def CMFtime_2_datetime(t):
    # internal CMF time is in minutes since 1850
    return datetime(1850, 1, 1) + timedelta(minutes = t)

# class DFM_model(_model):
#     def __init__(self, path_to_model, config_fn):
#         cmf_bmi = BMIWrapper(engine = path_to_model, config_fn = config_fn)
#         cmf_bmi.initialize()
#         super(CMF_model, self).__init__(cmf_bmi, 'Delft3D-FM', 'hydrodynamic')
#
#     def update(self, runoff, *args, **kwargs):
#         self.set_var('runoff', runoff)
#         self.bmi.update(*args, **kwargs)

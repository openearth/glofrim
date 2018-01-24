# -*- coding: utf-8 -*-
import numpy as np
from pcrglobwb_bmi_v203 import pcrglobwb_bmi
from bmi.wrapper import BMIWrapper
import logging

logger = logging.getLogger(__name__)

# wrapper around BMI
class _model(object):
    def __init__(self, bmi, name, mtype, missing_value=np.nan):
        self.bmi = bmi
        self.name = name
        self.mtype = mtype
        self.missing_value = missing_value
        # self.start_time = self.bmi.ge()

    def get_var(self, var, parse_missings=True):
        var_data = self.bmi.get_var(var)
        if parse_missings:
            # if given nodata is parsed to np.nan
            var_data = np.where(var_data == self.missing_value, np.nan, var_data)
        return var_data

    def set_var(self, *args, **kwargs):
        self.bmi.set_var(*args, **kwargs)

    def spinup(self, *args, **kwargs):
        self.bmi.spinup(*args, **kwargs)

    def update(self, *args, **kwargs):
        self.bmi.update(*args, **kwargs)
        current_time = self.bmi.get_current_time()
        time_step = self.bmi.get_time_step()
        logger.info(
            "start_time: %s, current_time %s, timestep %s",
            current_time, #TODO replace by self.start_time
            current_time,
            time_step
        )


class PCR_model(_model):
    def __init__(self, configfile, missing_value):
        pcr_bmi = pcrglobwb_bmi.pcrglobwbBMI()
        pcr_bmi.initialize(configfile)
        super(PCR_model, self).__init__(pcr_bmi, 'PCRGLOB-WB', 'hydrology',
                                        missing_value=missing_value)


class CMF_model(_model):
    def __init__(self, path_to_model, configfile, map_dir):
        cmf_bmi = BMIWrapper(engine = path_to_model, configfile = configfile)
        cmf_bmi.initialize(map_dir)
        super(CMF_model, self).__init__(cmf_bmi, 'CaMa-Flood', 'hydrodynamic')

    def update(self, runoff, *args, **kwargs):
        runoff = np.where(np.isnan(runoff), 0, runoff)
        self.set_var('runoff', runoff)
        super(CMF_model, self).update(*args, **kwargs)


# class DFM_model(_model):
#     def __init__(self, path_to_model, configfile):
#         cmf_bmi = BMIWrapper(engine = path_to_model, configfile = configfile)
#         cmf_bmi.initialize()
#         super(CMF_model, self).__init__(cmf_bmi, 'Delft3D-FM', 'hydrodynamic')
#
#     def update(self, runoff, *args, **kwargs):
#         self.set_var('runoff', runoff)
#         self.bmi.update(*args, **kwargs)

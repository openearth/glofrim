# -*- coding: utf-8 -*-

import numpy as np


def PCR2CMF_runoff(PCR_BMI, CMF_bmi):
    "creates a 'get_runoff' function to get runoff from PCR model fit for CMF model"
    def update_states():
        "coupling runoff between CMFan PCR model"
        runoff = PCR_BMI.get_var('landSurfaceRunoff')
        runoff = np.where(np.isnan(runoff), 0, runoff)
        CMF_bmi.set_var("runoff", runoff)
    return update_states

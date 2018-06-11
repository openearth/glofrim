# -*- coding: utf-8 -*-

import numpy as np

def get_CMF_delta_vol(CMF_bmi):
    # coupled 2d mask
    mask = (CMF_bmi.coupled_mask > 0).astype(float)
    # total Q inflow 2D for coupled cells
    rivinf = np.copy(CMF_bmi.get_var('rivinf'))
    rivinf = np.where(np.isnan(rivinf), 0, rivinf)
    fldinf = np.copy(CMF_bmi.get_var('fldinf'))
    fldinf = np.where(np.isnan(fldinf), 0, fldinf)
    q_in = (rivinf + fldinf) * mask # [m3/s]
    # total CMF converted runoff inflow 2D for coupled cells
    runoff = np.copy(CMF_bmi.get_var('runoff'))
    runoff = np.where(np.isnan(runoff), 0, runoff)
    runoff = runoff * mask # [m3/s]
    # take Qin + runoff in most upstream coupled cells, only runoff for other cells
    tot_flux = np.where(CMF_bmi.coupled_mask == 2, q_in + runoff, runoff)
    # convert flux to volume per day
    delta_vol = tot_flux * CMF_bmi.options['dt'] # [m3/day]
    return delta_vol

def set_CMF_forcing(PCR_bmi, CMF_bmi):
    """coupling runoff between CMF and PCR model"""
    runoff = np.copy(PCR_bmi.get_var('runoff')) # [m/day]
    runoff = np.where(np.isnan(runoff), 0, runoff)
    # note that runoff in (roffin) should be used to set PCR runoff.
    CMF_bmi.set_var("roffin", runoff) # [m/dtin] = [m/d]

def set_DFM_forcing(DFM_bmi, CMF_bmi):
    # determine coupled indices, cell area in DFM
    DFMidx = DFM_bmi.coupled_idx
    DFM_area_1d = DFM_bmi.get_var('ba')[DFMidx]
    # determine coupled indices and cell area fraction in CMF
    CMFidx = CMF_bmi.coupled_idx
    CMFfrac = CMF_bmi.coupled_area_frac
    # get CMF delta volume
    CMF_delta_vol = get_CMF_delta_vol(CMF_bmi) # 2d array [m3/day]
    # convert to fit to DFM rain variable unitss
    DFM_depth_conservative = CMF_delta_vol[CMFidx] * CMFfrac / DFM_area_1d # 1d array for coupld DFM cells [m/day]
    zerorain = np.zeros_like(np.copy(DFM_bmi.get_var('rain')))
    DFM_bmi.set_var('rain', zerorain)
    DFM_bmi.set_var_index('rain', DFMidx, DFM_depth_conservative)

# def PCR2CMF_runoff(PCR_BMI, CMF_bmi):
#     "creates a 'get_runoff' function to get runoff from PCR model fit for CMF model"
#     def update_states():
#         "coupling runoff between CMFan PCR model"
#         runoff = PCR_BMI.get_var('landSurfaceRunoff')
#         runoff = np.where(np.isnan(runoff), 0, runoff)
#         CMF_bmi.set_var("runoff", runoff)
#     return update_states

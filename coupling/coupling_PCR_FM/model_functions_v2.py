# -*- coding: utf-8 -*-
import numpy as np
from pcrglobwb_bmi_v203 import pcrglobwb_bmi
from bmi.wrapper import BMIWrapper
import logging
import warnings
from datetime import datetime, timedelta
import os, shutil
logger = logging.getLogger(__name__)

# local libraries
from utils import write_ini, set_values_in_array
from coupling_functions import assignPCR2cells

# wrapper around BMI
class _model(object):
    def __init__(self, bmi, name, t_unit, missing_value=np.nan):
        self.bmi = bmi
        self.name = name
        self.t_unit = t_unit
        self.missing_value = missing_value
        self.dt=-1 # for initialization only. must be overwritten by initialize function later

    # model initialization
    def initialize(self, *args, **kwargs):
        self.bmi.initialize(*args, **kwargs)
        # set start time attribute
        self.start_time = self.get_start_time()

    def spinup(self, *args, **kwargs):
        self.bmi.spinup(*args, **kwargs)

    # set model and coupling options
    def set_options(self, dt=None, *arg, **kwargs):
        if dt is not None:
            self.dt = dt
            logger.info('{:s} model dt set to {:.2f} {:s}'.format(
                                        self.name, self.dt, self.t_unit))

    def set_update_states(self, update_states):
        if not callable(update_states):
            raise Warning("update_states argument should be a callable function")
        self.update_states = update_states

    # exchange states
    def get_var(self, var, parse_missings=True, *args, **kwargs):
        var_data = self.bmi.get_var(var)
        if parse_missings:
            # if given nodata is parsed to np.nan
            var_data = np.where(var_data == self.missing_value, np.nan, var_data)
        return var_data

    def set_var(self, *args, **kwargs):
        self.bmi.set_var(*args, **kwargs)

    # update states placeholder
    def update_states(self):
        warning.warn("no update states function set",
                     RuntimeWarning
                     )

    # run timestep dt in model
    def update(self, dt=None):
        if dt is None: # by default take internally set dt
            dt = self.dt
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

    # time attributes
    def get_start_time(self):
        "get model start (initialization) time"
        return self.bmi.get_start_time()

    def get_current_time(self):
        "get model current time"
        return self.bmi.get_current_time()

    def get_time_step(self):
        "get model current time step"
        return self.bmi.get_time_step()



class PCR_model(_model):
    def __init__(self, missing_value=-999, landmask_mv=255):
        # BMIWrapper for PCR model model
        pcr_bmi = pcrglobwb_bmi.pcrglobwbBMI()
        super(PCR_model, self).__init__(pcr_bmi, 'PCRGLOB-WB', 'day',
                                        missing_value=missing_value)
        self.landmask_mv = landmask_mv

    def initialize_offline_forcing(self, config_fn,
                                  start_date, end_date, out_dir, in_dir,
                                  dt=1, ini_kwargs={}):
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        self.out_dir = out_dir

        # make ini file consistent with GLOFRIM settings
        date_fmt = "%Y-%m-%d"
        ini_kwargs.update({
            "inputDir": os.path.abspath(in_dir),
            "outputDir": os.path.abspath(out_dir),
            "startTime": start_date.strftime(date_fmt),
            "endTime": end_date.strftime(date_fmt),
            })
        tmp_config_fn = os.path.join(out_dir, 'tmp_' + os.path.basename(config_fn))
        write_ini(tmp_config_fn, config_fn, **ini_kwargs)
        logger.info('Ini file for PCR model written to {:s}'.format(tmp_config_fn))
        self.config_fn = tmp_config_fn

        # NOTE: this changes cwd and might cause errors downstream
        super(PCR_model, self).initialize(config_file_location=tmp_config_fn)
        logger.info('PCR model initialized')
        self.set_options(dt=dt)

    def get_model_coords(self):
        # TODO read transform
        transform = None
        self.index = index
        #self.center_points_2d =

    def get_total_water_volume(self): #model_pcr, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr):
        """
        Calculating the delta volumes [m3/d] for PCR-cells.
        Delta volumes are based on discharge, surfaceRunoff, and topWaterLayer (only 2way-coupling) in PCR-GLOBWB.

        Input:
            - list with indexes pointing to coupled PCR-cells
            - PCR landmask data

        Output:
            - if RFS active, two arrays with delta volumes for river and floodplain cells, respectively
            - if RFS not active, one array with aggregated delta volumes
            - all outputs are in m3/day
        """
        #- retrieve data from PCR-GLOBWB
        current_discharge_pcr  = self.get_var('discharge')
        current_runoff_pcr     = self.get_var('landSurfaceRunoff')
        current_waterlayer_pcr = self.get_var('topWaterLayer')

        # 1a. Discharge
        # loop over current discharge and convert to m3/d; missing values are replaced with zero
        water_volume_PCR_rivers = current_discharge_pcr * 86400. #sec/day
        water_volume_PCR_runoff = current_runoff_pcr * self.get_var('cellArea')
        water_volume_PCR_waterlayer = current_waterlayer_pcr * self.get_var('cellArea')
        total_water_volume = water_volume_PCR_rivers + water_volume_PCR_runoff + water_volume_PCR_waterlayer

        return total_water_volume

    def deactivate_LDD(self, indices):
        """
        NOTE:
        :param indices: list of (y, x) indices / 'all'
        """
        # retrieving current LDD map
        LDD_PCR_new = np.copy(self.get_var(('routing', 'lddMap'), parse_missings=False))
        if indices != 'all':
            # replace LDD values within the hydrodynamic model domain to 5 (pit cell)
            LDD_PCR_new = set_values_in_array(LDD_PCR_new, indices, 5)
        else:
            LDD_PCR_new = np.where(LDD_PCR_new != self.landmask_mv, 5, self.landmask_mv)
        # overwriting current with new LDD information
        self.set_var(('routing', 'lddMap'), LDD_PCR_new, self.landmask_mv)



class CMF_model(_model):
    def __init__(self, engine,
                 missing_value=1e20):
        """GLOFRIM wrapper around BMI for CaMa-Flood (CMF)

        input
        engine          BMIWrapper path to CMF engine
        config_fn      BMIWrapper CMF config file name
        map_dir         CMF map directory. Required to initialize CMF
        update_states   Function to get update the states (e.g. runoff) based on coupled model.
                        This can also be set later using the set_options function.

        """
        ## initialize BMIWrapper and model
        # NOTE that that the CMF BMI wrapper does not use the configfile argument
        cmf_bmi = BMIWrapper(engine = engine)
        super(CMF_model, self).__init__(cmf_bmi, 'CaMa-Flood', 'sec',
                                        missing_value=missing_value)

    def initialize_online_forcing(self, config_fn, model_dir,
                                  start_date, end_date, out_dir,
                                  dt=86400, ini_kwargs={}):

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        # write ini file
        ini_kwargs.update(self._ini_kwargs(out_dir, start_date, end_date))
        ini_kwargs.update({"CRUNOFFDIR": "",
                           "LBMIROF": ".TRUE.",
                           "DTIN": "86400", # input unit settings -> runoff in m/s
                           "DROFUNIT": "1",
                           })
        # TODO: write consistent inpmat file and add args below to ini_kwargs
        # "CDIMINFO":  # input dimensions file
        # "CINPMAT":

        # TODO: write input_flood.nam to out_dir and adapt model dir.
        # TODO: add sym link in out_dir to map dir to keep relative paths
        tmp_config_fn = os.path.join(model_dir, "input_flood.nam")
        if os.path.abspath(tmp_config_fn) == os.path.abspath(config_fn):
            shutil.move(config_fn, config_fn + '.template')
            config_fn = config_fn + '.template'
            logger.warning('template .nam file cannot be called input_flood.nam; file renamed {:s}'.format(config_fn))
        write_ini(tmp_config_fn, config_fn, ignore='!', **ini_kwargs)
        self.config_fn = tmp_config_fn
        logger.info('tmp ini file for CMF model written to {:s}'.format(tmp_config_fn))

        # initialize model
        self.initialize(configfile=model_dir)
        logger.info('CMF model initialized')
        self.set_options(dt=dt)

    def get_model_coords(self):
        # TODO read transform
        transform = None
        self.transform = transform

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

    # internal function
    def _ini_kwargs(self, out_dir, start_date, end_date, **kwargs):
        ## input_flood.nam key word arguments that alwasys need to be adapted
        ini_kwargs = {"COUTDIR": os.path.abspath(out_dir), # output dir
                      "ISYEAR": "{:d}".format(start_date.year),  # model start & end dates
                      "ISMON": "{:d}".format(start_date.month),
                      "ISDAY": "{:d}".format(start_date.day),
                      "IEYEAR": "{:d}".format(end_date.year),
                      "IEMON": "{:d}".format(end_date.month),
                      "IEDAY": "{:d}".format(end_date.day),
                      }
        return ini_kwargs


class DFM_model(_model):
    def __init__(self, engine):
        dfm_bmi = BMIWrapper(engine = engine)
        super(DFM_model, self).__init__(dfm_bmi, 'Delft3D-FM', 'sec')

    def initialize_online_forcing(self, config_fn, model_dir,
                                  start_date, end_date, out_dir,
                                  dt=86400, ini_kwargs={}):

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        # write ini file
        # TODO: make required changes to .mdu file
        ini_kwargs.update({})
        tmp_config_fn = os.path.join(model_dir, 'tmp_' + os.path.basename(config_fn))
        write_ini(tmp_config_fn, config_fn, **ini_kwargs)
        self.config_fn = tmp_config_fn
        logger.info('tmp ini file for DFM model written to {:s}'.format(tmp_config_fn))

        # initialize model
        self.initialize(configfile=tmp_config_fn)
        logger.info('CMF model initialized')
        self.set_options(dt=dt)

        # get Model coordinates
        self.get_model_coords()

    def get_model_coords(self):
        # define separator between 2D and 1D parts of arrays == lenght of 2d cell points
        self.separator_1d2d = len(self.get_var('flowelemnode'))
        x_coords = self.get_var('xz') # x-coords of each cell centre point
        y_coords = self.get_var('yz') # y-coords of each cell centre point
        xy_coords = [[(xy)] for xy in zip(x_coords, y_coords)]
        self.center_points_2d = xy_coords[:self.separator_1d2d]
        self.center_points_1d = xy_coords[self.separator_1d2d:]


    def calculateDeltaWater(self, total_water_volume, af_list):
        delta_water = []
        for i in xrange(len(total_water_volume)):
            added_water_level = total_water_volume[i] / af_list[i]
            delta_water.extend(added_water_level)

    def set_var_1d(self, var, data):
        zeroValuesFor2D  = np.zeros(self.separator_1d2d)
        super(DFM_model, self).set_var(var, np.stack([zeroValuesFor2D, data]))

    # def get_var_1d(var, *args, **kwargs):
    #     var = super(DFM_model, self).get_var(var, *args, **kwargs)
    #     return var[self.separator_1d2d:]
    #
    # def get_var_2d(var, *args, **kwargs):
    #     var = super(DFM_model, self).get_var(var, *args, **kwargs)
    #     return var[:self.separator_1d2d]


# utils
def CMFtime_2_datetime(t):
    # internal CMF time is in minutes since 1850
    return datetime(1850, 1, 1) + timedelta(minutes = t)

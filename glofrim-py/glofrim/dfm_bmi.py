import numpy as np
import sys
from configparser import ConfigParser
import logging
import os
from os.path import join, isfile, abspath, dirname, basename, relpath
from datetime import datetime, timedelta
import rasterio
import re

from bmi.wrapper import BMIWrapper as _bmi

from utils import setlogger
from gbmi import GBmi
from grids import UGrid
import glofrim_lib as glib 

class DFM(GBmi):
    """
    Glofrim implementation of the PCR BMI adaptor.
    """
    _name = 'DFM'
    _long_name = 'Delf3D-FM'
    _version = 'x.x.x'
    _var_units = {'rain': 'mm', 'ba': 'm2'}
    _input_var_names = ['rain', 'ba']
    _output_var_names = ['hs', 's1']
    _area_var_name = 'ba'
    _tdict = {'H': 'hours', 'M': 'minutes', 'S': 'seconds'} # mdu time unit to datatime format

    def __init__(self, engine):
        self._bmi = _bmi(engine = engine)
        self.logger = setlogger(None, self._name, thelevel=logging.INFO)
        self.initialized = False
        self.grid = None

    """
    Model Control Functions
    """
    def initialize_config(self, config_fn):
        if self.initialized:
            raise Warning("model already initialized, it's therefore not longer possible to initialize the config")
        # config settings
        self._config_fn = abspath(config_fn)
        self._config = glib.configread(self._config_fn, encoding='utf-8', cf=ConfigParser(inline_comment_prefixes=('#')))
        self._datefmt = "%Y%m%d"
        # model time
        self._timeunit = self.get_time_units()
        self._dt = self.get_time_step()
        self._startTime = self.get_start_time()
        self._endTime = self.get_end_time()
        self._t = self._startTime
        # model files
        _root = dirname(self._config_fn)
        self._outdir = glib.getabspath(self.get_attribute_value('output:OutputDir'), _root)
        self.logger.info('Config initialized')

    def initialize_model(self, **kwargs):
        if not hasattr(self, '_config_fn'):
            raise Warning('Run initialize_config before initialize_model')
        self.write_config() # write updated config to file as bmi does not allow direct access
        self._bmi.initialize(self._config_fn)
        self.initialized = True
        self.logger.info('Model initialized')

    def initialize(self, config_fn):
        if not hasattr(self, '_config'):
            self.initialize_config(config_fn)
        self.initialize_model()
            
    def update(self):
        if self._t >= self._endTime:
		    raise Exception("endTime already reached, model not updated")
        self._bmi.update(dt=self.get_time_step().days)
        self._t += self.get_time_step()
        self.logger.info('updated model to datetime {}'.format(self._t.strftime("%Y-%m-%d")))

    def update_until(self, t):
        if (t<self._t) or t>self._endTime:
            raise Exception("wrong time input: smaller than model time or larger than endTime")
        while self._t < t:
            self.update()

    # not defined in CMF
    def spinup(self):
        """PCR specific spinup function"""
        raise NotImplementedError()

    def finalize(self):
        self._bmi.finalize()


    """
    Variable Information Functions
    """

    def get_start_time(self):
        RefDate = self.get_attribute_value('time:RefDate')
        RefDate = datetime.strptime(RefDate, self._datefmt)
        if self.initialized:
            # date to datetime object
            TStart = self._bmi.get_start_time()
        else:
            TStart = float(self.get_attribute_value('time:TStart'))
        startTime = RefDate + timedelta(**{self.get_time_units(): TStart})
        self._startTime = startTime
        return self._startTime
    
    def get_current_time(self):
        if self.initialized:
            return self._bmi.get_current_time()
        else:
            return self.get_start_time()

    def get_end_time(self):
        RefDate = self.get_attribute_value('time:RefDate')
        RefDate = datetime.strptime(RefDate, self._datefmt)
        if self.initialized:
            # date to datetime object
            TStop = self._bmi.get_end_time()
        else:
            TStop= float(self.get_attribute_value('time:TStop'))
        endTime = RefDate + timedelta(**{self.get_time_units(): TStop})
        self._endTime = endTime
        return self._endTime

    def get_time_step(self):
        if not hasattr(self, '_dt'):
            dt = float(self.get_attribute_value('time:DtUser'))
            self._dt = timedelta(**{self.get_time_units(): dt})
        return self._dt 
        
    def get_time_units(self):
        if not hasattr(self, '_timeunit'):
            tunit = self.get_attribute_value('time:Tunit')
            if not tunit in self._tdict.keys():
                msg = 'Unknown time unit in time:Tunit attribute: {}. Use {}.'
                raise ValueError(msg.format(tunit, ', '.join(self._tdict.keys())))
            self._timeunit = self._tdict[tunit]
        return self._timeunit


    """
    Variable Getter and Setter Functions
    """
    
    def get_value(self, long_var_name, **kwargs):
        return np.asarray(self._bmi.get_var(long_var_name))

    def get_value_at_indices(self, long_var_name, inds, **kwargs):
        return self.get_value(long_var_name).flat[inds]

    def set_value(self, long_var_name, src, **kwargs):
        self._bmi.set_var(long_var_name, src)

    def set_value_at_indices(self, long_var_name, inds, src, **kwargs):
        # always use 1d inds
        if not ((isinstance(inds, np.ndarray)) and (inds.ndim == 1)):
            raise ValueError('indices should be 1d arrays')
        val = self.get_value(long_var_name)
        val.flat[inds] = src
        self._bmi.set_var(long_var_name, val)

    """
    Grid Information Functions
    """

    def get_grid(self):
        """Get DFM model coordinates for 1D and 2D mesh via BMI. The DFM model
        should be initialized first in order to access the variables."""
        if not hasattr(self, 'grid') or (self.grid is None):
            self.logger.info('Getting DFM UGrid.')
            # all cells 2d + 1d
            cell_x = self.get_value('xz') # x-coords of each cell centre point
            cell_y = self.get_value('yz') # y-coords of each cell centre point
            cell_xy = np.array(zip(cell_x, cell_y))
            # 1D-2D divide for cells
            n2d = int(self.get_value('ndx2d'))
            cidx_1d = np.arange(n2d, cell_xy.shape[0])
            # all nodes first 2D, then 1D, then boundaries
            node_lon = self.get_value('xk')
            node_lat = self.get_value('yk') 
            nodes = np.array(zip(node_lon, node_lat))
            # 2d mesh 
            faces = self.get_value('flowelemnode') # face_node_connectivity
            nidx_2d = np.arange(faces.max()) # index of 2d nodes
            faces = np.ma.masked_equal(faces, -1) - 1
            nodes_2d = nodes[nidx_2d, :]     # coordinates of nodes
            face_coordinates = cell_xy[:n2d, :] # coordinates of cell centres
            self.grid = UGrid(nodes_2d, faces, face_coordinates=face_coordinates, fill_value=-1, crs=None)
            ## 1D network
            # TODO: indices of nodes and cell-centers (at which data is stored) do not match for 1d ?? check with Arthur v. Dam
            # if we figure out how the cell indices relate to the node indices we can set a 1d network with links.
            self.grid.set_1d(cell_xy[cidx_1d, :], links=None)
        return self.grid


    # def get_snetwork(self):
    #     """Get DFM 2D mesh definition via BMI"""
    #     if not hasattr(self, 'snetwork') or (self.snetwork is None):
    #         import rtree
    #         self.logger.info('Getting DFM 1D network.')
    #         # all cells 2d + 1d
    #         cell_x = self.get_value('xz') # x-coords of each cell centre point
    #         cell_y = self.get_value('yz') # y-coords of each cell centre point
    #         cell_xy = np.array(zip(cell_x, cell_y))
    #         # 1D-2D divide for cells
    #         n2d = int(self.get_value('ndx2d'))
    #         cidx_1d = np.arange(n2d, cell_xy.shape[0])
    #         # all nodes first 2D, then 1D, then boundaries
    #         node_lon = self.get_value('xk')
    #         node_lat = self.get_value('yk') 
    #         nodes = np.array(zip(node_lon, node_lat))
    #         # all links
    #         kn = self.get_value('kn') 
    #         nlink_type = kn[:, 2] # type
    #         nlinks = kn[:, :2] - 1 # link
    #         nlinks_1d = nlinks[nlink_type == 1] # link_between_1D_nodes
    #         nidx_1d = np.arange(nlinks_1d.min(), nlinks_1d.max()+1)
    #         nodes_1d = nodes[nidx_1d, :]
    #         nlinks_1d = nlinks_1d - nlinks_1d.min()
    #         # create index
    #         ncidx_1d = np.ones_like(nidx_1d) * -1
    #         rt = rtree.index.Index()
    #         for ci, xy in zip(cidx_1d, cell_xy[cidx_1d]):
    #             rt.insert(ci, xy.tolist())
    #         # find nearest cell for each node
    #         for ni, xy in enumerate(nodes[nidx_1d]):
    #             ci = list(rt.nearest(xy.tolist(), 1))
    #             assert np.all(xy == cell_xy[ci])
    #             ncidx_1d[ni] = ci[0] 
    #             # make sure each cell is found only ones
    #             # this is expensive, but only way to do it
    #             rt.delete(ci[0], xy.tolist()) 
    #         self.snetwork = SNetwork(nodes_1d, nlinks_1d, inds=ncidx_1d, fill_value=-1, crs=None)
    #     return self.snetwork

    """
    set and get attribute / config 
    """

    def set_start_time(self, start_time):
        if isinstance(start_time, datetime):
            RefDate = start_time.strftime(self._datefmt)
            TStart = timedelta(hours=start_time.hour, minutes=start_time.minute, seconds=start_time.second).seconds
            if self.get_time_units() == 'minutes':
                TStart = TStart / 60
            elif self.get_time_units() == 'hours':
                TStart = TStart / 3600
            TStart = '{:.0f}'.format(TStart)
        else:
            TStart = '0.'
        try:
            start_time = datetime.strptime(start_time, self._datefmt) 
            self._startTime = start_time
        except ValueError:
            raise ValueError('wrong date format, use "yyyymmdd"')
        self.set_attribute_value('time:RefDate', RefDate)
        self.set_attribute_value('time:TStart', TStart)

    def set_end_time(self, end_time):
        if isinstance(end_time, str):
            try:
                end_time = datetime.strptime(end_time, self._datefmt)
            except ValueError:
                raise ValueError('wrong end_date format, use "yyyymmdd"')
        if not isinstance(end_time, datetime):
            raise ValueError('wrong end_date datatype')
        RefDate = self.get_attribute_value('time:RefDate')
        RefDate = datetime.strptime(RefDate, self._datefmt)
        assert end_time >  RefDate
        self._endTime = end_time
        TStop = (end_time - RefDate).seconds + (end_time - RefDate).days * 86400
        if self.get_time_units() == 'minutes':
            TStop = TStop / 60
        elif self.get_time_units() == 'hours':
            TStop = TStop / 3600
        TStop = '{:.0f}'.format(TStop)
        self.set_attribute_value('time:TStop', TStop)

    def set_out_dir(self, out_dir):
        self.set_attribute_value('output:OutputDir', abspath(out_dir))
        self._outdir = abspath(out_dir)

    def get_attribute_names(self):
        glib.configcheck(self, self.logger)
        return glib.configattr(self._config)
    
    def get_attribute_value(self, attribute_name):
        glib.configcheck(self, self.logger)
        self.logger.debug("get_attribute_value: {}".format(attribute_name))
        return glib.configget(self._config, attribute_name)
    
    def set_attribute_value(self, attribute_name, attribute_value):
        glib.configcheck(self, self.logger)
        self.logger.debug("set_attribute_value: {} -> {}".format(attribute_name, attribute_value))
        return glib.configset(self._config, attribute_name, str(attribute_value))

    def write_config(self):
        """write adapted config to file. just before initializing
        only for models which do not allow for direct access to model config via bmi"""
        glib.write_config(self, self._config, self._config_fn, self.logger)
from glofrim.grids import GridType
import os
from os.path import isfile, relpath, dirname
import glofrim.glofrim_lib as glib
import rasterio.warp
import rasterio.enums

import numpy as np

class SpatialCoupling(object):
    _methods = ['grid_2_grid', 'grid_2_1d', 'grid_us_2_1d_us']

    def __init__(self, to_bmi, from_bmi, method=None, filename=None, to_coords=None):
        # this info can be filled now or when coupling
        self.method = method
        if not self.method in self._methods:
            raise AttributeError('Coupling method {} does not exist'.format(method))
        self.filename = filename
        if self.filename is not None and not isfile(self.filename):
            raise IOError('{} not found.'.format(self.filename))
        # from
        self.from_bmi = from_bmi
        self.from_name = from_bmi.get_component_name()
        # to
        self.to_bmi = to_bmi
        self.reproject = None
        self.to_name = to_bmi.get_component_name()
        # this info is filled when coupling
        self.from_ind = np.array([]) # 1d array with indices
        self.from_grp = np.array([]) # 1d array with group number
        self.from_grp_n = np.array([]) # 1d array with indices per group
        self.to_ind = np.array([])
        self.to_grp = np.array([])
        self.to_grp_n = np.array([])
        self.to_coords = to_coords
        self.frac = np.array([])
        self._nodes_outside_domain = []

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        s = "SpatialCoupling object from {} to {} based on {} method"
        return s.format(self.from_name, self.to_name, self.method)

    def set_inds_from_coupled_dict(self, coupled_dict):
        if -1 in coupled_dict.keys():
            not_coupled = coupled_dict.pop(-1)
            print(not_coupled)
            # TODO warning
        self.to_ind, self.to_grp, self.to_grp_n = group_index(coupled_dict.values())
        self.from_ind, self.from_grp, self.from_grp_n = group_index(coupled_dict.keys()) 
    
    def set_frac(self, from_area=True):
        if not from_area:
            self.frac = np.ones_like(self.to_ind).astype(float) / np.repeat(self.to_grp_n, self.to_grp_n)
        else:
            area_var = self.to_bmi._area_var_name
            areas = self.to_bmi.get_value_at_indices(area_var, self.to_ind)
            area_sum = groupby_sum(areas, self.to_grp) 
            self.frac = areas.astype(float) / np.repeat(area_sum, self.to_grp_n)
        assert np.all(groupby_sum(self.frac, self.to_grp).round(5)==1)
        return self.frac
            
    def get_frac(self, from_area=True):
        if (self.to_ind.size > 0) and (self.frac.size == 0):
            self.set_frac(from_area=from_area)
        return self.frac

    def couple(self, method=None, filename=None, **kwargs):
        # replace properties if provided (not None)
        self.method = method if method is not None else self.method
        if not self.method in self._methods:
            raise AttributeError('Coupling method {} does not exist'.format(method))
        self.filename = filename if filename is not None else self.filename
        if self.filename is not None and not isfile(self.filename):
            raise IOError('{} not found.'.format(self.filename))
        # do coupling
        self.from_bmi.get_grid() # make sure grid property is set
        self.to_bmi.get_grid() # make sure grid property is set
        getattr(self, self.method)()

    def to_file(self, filename):
        # TODO: save to JSON
        raise NotImplementedError()

    def from_file(self, filename):
        # TODO: read from JSON
        raise NotImplementedError()

    def grid_2_grid(self):
        """
        """
        # rgrid to rgrid
        if (self.to_bmi.grid.type == 1) and (self.from_bmi.grid.type == 1):
            # just check grid types. do nothing
            # TODO make a configurable Resampling choice, now Resampling.average used hard coded
            check_bounds= np.all(self.to_bmi.grid.bounds == self.from_bmi.grid.bounds)
            check_shape = np.all(self.to_bmi.grid.shape == self.from_bmi.grid.shape)
            if not (check_shape and check_bounds):
                # make a reprojection function, note that rasterio sees rasters flipped 
                # vertically with respect to BMI convention, so a flipud is required
                self.reproject = lambda data, nodata: rasterio.warp.reproject(
                    data,
                    destination=np.zeros((self.to_bmi.grid.height, self.to_bmi.grid.width)),
                    src_transform=self.from_bmi.grid.transform,
                    src_crs=self.from_bmi.grid.crs,
                    src_nodata=nodata,
                    dst_transform=self.to_bmi.grid.transform,
                    dst_crs=self.to_bmi.grid.crs,
                    dst_nodata=nodata,
                    resampling=rasterio.enums.Resampling.average
                    )[0]
                # raise ValueError('both model grids should have the shape and bounding box')
            # TODO make special case to cover same rgrid and same ugrid coupling
        # rgrid to ucgrid
        elif (self.from_bmi.grid.type == 1) and (self.to_bmi.grid.type == 3):
            # set inpmat
            bounds, res = self.from_bmi.grid.bounds, self.from_bmi.grid.res
            olat = 'NtoS' if self.from_bmi.grid.NtoS else 'StoN'
            self.to_bmi.set_inpmat(bounds, res, olat=olat)
            # set new inpmat and diminfo in config
            self.to_bmi.set_inpmat_attrs()
        else:
            # TODO add options to couple u/rgrid to u/rgrid based on cell (face) centers
            # possibilities: 
            # 1) to_grid < from_grid: find from_cell for each to_cell. check from cells which are not touched
            # 2) to_grid > from_grid: raise warning?? same procedure will do. 

            msg = 'coupling from grid type {} to grid type {} has not been implemented'
            raise NotImplementedError(msg.format(self.to_bmi.grid.type, self.from_bmi.grid.type))

    def grid_2_1d(self):
        """
        """
        # no working index yet for UGrid
        if self.from_bmi.grid.type == 2:
            msg = 'coupling from grid type {} to 1d has not been implemented'
            raise NotImplementedError(msg.format(self.from_bmi.grid.type))
        if not self.to_bmi.grid.has_1d():
            raise ValueError('1d network not set for model {}'.format(self.to_name))
        # get 1d nodes
        to_coords = self.to_bmi.grid._1d.nodes
        to_inds = self.to_bmi.grid._1d.inds
        # find cells indices for all nodes
        from_inds = self.from_bmi.grid.index(to_coords[:, 0], to_coords[:, 1])
        valid = from_inds >= 0
        if np.any(~valid):
            print('1D nodes found outside of valid 2D domain')
            self._nodes_outside_domain = to_inds[~valid]
        if np.all(~valid):
            raise IndexError('All 1D nodes outside of valid 2D domain')
        # to(1):from(1) dict
        d = dict(zip(to_inds[valid].tolist(), from_inds[valid].tolist()))
        # from(n):from(1) dict
        coupled_dict = dictinvert(d)
        # set indices
        self.set_inds_from_coupled_dict(coupled_dict)
        return coupled_dict


    def grid_us_2_1d_us(self):
        """"""
        if not self.from_bmi.grid.has_dd():
            raise ValueError('dd network not set for model  {}'.format(self.from_name))
        if not self.to_bmi.grid.has_1d():
            raise ValueError('1d network not set for model  {}'.format(self.to_name))
        # get 1d nodes
        if self.to_coords is None:
            to_coords = self.to_bmi.grid._1d.nodes
            to_inds = self.to_bmi.grid._1d.inds
        else:
            # get 1d nodes from user-configured list of coordinates in .ini file
            to_coords = np.array(self.to_coords)
            x, y = zip(*to_coords)
            row, col = rasterio.transform.rowcol(self.to_bmi.grid.transform, x, y)
            try:
                to_inds = np.array(self.to_bmi.grid.ravel_multi_index(row, col))
            except:
                raise IndexError('Some or all of manually configured 1D nodes fall outside valid model domain')

        # get short handles for function
        fup = self.from_bmi.grid._dd.next_upstream
        findex = self.from_bmi.grid.ravel_multi_index
        frowcol = self.from_bmi.grid.unravel_index
        # find cells indices for all nodes
        from_inds = self.from_bmi.grid.index(to_coords[:, 0], to_coords[:, 1], src_crs=self.to_bmi.grid.crs)
        valid = from_inds >= 0
        if np.any(~valid):
            print('1D nodes found outside of valid 2D domain')
            self._nodes_outside_domain = to_inds[~valid]
        if np.all(~valid):
            raise IndexError('All 1D nodes outside of valid 2D domain')
        # find  unique indices and convert to row, col
        u_from_inds = np.unique(from_inds[valid])
        from_rows, from_cols = frowcol(from_inds[valid])
        # to(1):from(n) dict
        d = dict()
        for r, c, i in zip(from_rows, from_cols, to_inds[valid]):
            from_idx_us = findex(*fup(r, c)) # upstream cell of node i
            # find most upstream of connected cells
            from_idx_us = [idx for idx in from_idx_us if idx not in u_from_inds]
            if len(from_idx_us) > 0:
                d[i] = tuple(from_idx_us) # tuple as list cannot be used as python dict key
        # from(n):from(n) dict
        coupled_dict = dictinvert(d)
        if len(d) == 0:
            raise IndexError('There are no 2D cells upstream from most upstream 1D nodes')
        # set indices
        self.set_inds_from_coupled_dict(coupled_dict)
        return coupled_dict

    def at_indices(self):
        """check if coupling at indices"""
        return self.from_ind.size > 0

def dictinvert(d):
    inv = {}
    for k, v in d.items():
        keys = inv.setdefault(v, [])
        keys.append(k)
    return inv

def group_index(indices):
    """
    translate jagged array with indices to flattened array and groups
    
    group_index([[0,1,2], 
                 [3,4]])
    returns
    index:  [0, 1, 2, 3, 4]
    groups: [0, 0, 0, 1, 1]
    length: [3, 2]
    """
    x = np.array(list(indices)).flatten()
    if isinstance(x[0], (list, tuple)): # check if jagged array
        grp_n = np.array([len(a) for a in x])
        grp = np.repeat(np.arange(x.size), grp_n)
        x = np.concatenate(x)
    else:
        grp_n = np.ones(x.size)
        grp = np.arange(x.size, dtype=np.int32)
    return x.astype(np.int32), grp, grp_n

def groupby_sum(X, groups):
    """vectorized sum of groups in X"""
    uf = np.add
    X = np.array(X).flatten()
    groups = np.array(groups).flatten()
    if not X.size == groups.size:
        raise ValueError('array and group sizes do not match')
    minlength = groups.max() + 1
    out = np.full((minlength,), uf.identity, dtype = X.dtype)
    uf.at(out, groups, X)
    return out

# def jagged2masked_idx(data, fill_value=-1):
#     # Get lengths of each row of data
#     data = np.array(data)
#     lens = np.array([len(np.atleast_1d(i)) for i in data])
#     if np.all(lens == 1):
#         out = data[:, None]
#     else:
#         # Mask of valid places in each row
#         mask = np.arange(lens.max()) < lens[:,None]
#         # Setup output array and put elements from data into masked positions
#         out = np.ones(mask.shape, dtype=np.int32) * fill_value
#         out[mask] = np.concatenate(data)
#     return np.ma.masked_equal(out, fill_value)
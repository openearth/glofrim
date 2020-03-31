#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Dirk Eilander (contact: dirk.eilancer@vu.nl)
# Created: Nov-2017

import numpy as np
from .nb_ops import MaskedWdwArray
# TODO: add io, check if sound and bbox-clip options
# TODO: test if drainge area etc function work fast compared to pcraster


class DrainageDirection(object):
    """Drainage direction numpy-based object to do fast operations on (nextx, nexty),
    LDD or D8 drainge direction networks"""

    def __init__(self, stype, dd, raster, transform,
                 pit_value=0, fill_value=255,
                 search_radius=1, **kwargs):
        # general properties
        self.type = stype
        self.dd = dd
        self.pit_value = pit_value
        # create up and downstream nb patterns
        if dd is not None:
            self._dd_ds = self.dd.flatten()
            self._dd_us = np.fliplr(np.flipud(self.dd)).flatten()
            self._dd_ds = np.where(
                self._dd_ds == pit_value, pit_value-5, self._dd_ds)
            self._dd_us = np.where(
                self._dd_us == pit_value, pit_value-5, self._dd_us)
        else:  # in case of nextxy
            self._dd_ds, self._dd_us = None, None
        # initialize masked array
        raster = np.ma.masked_equal(raster, fill_value, copy=False)
        self.r = MaskedWdwArray(raster, sr=search_radius)
        self.shape = self.r.shape
        if self.r.ndim == 3:
            self._mask2d = ~self.r.mask[0, :, :]
        else:
            self._mask2d = ~self.r.mask
        # set geo properties
        self.transform = transform

    # general relations
    def _nb_downstream_idx(self, row, col):
        """neighborhood downstream index"""
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        valid = np.logical_not(self._pit(row, col), self.r.mask[row, col])
        row, col = row[valid], col[valid]
        wdw_bool = self._dd_ds[None, :] == self.r.data[row, col][:, None]
        wdw_rows, wdw_cols = self.r.wdw(row, col)[1]
        row_ds, col_ds = np.atleast_1d(
            wdw_rows[wdw_bool]), np.atleast_1d(wdw_cols[wdw_bool])
        return (row_ds, col_ds), valid

    def _nb_upstream_idx(self, row, col):
        """neighborhood upstream index"""
        (row_us, col_us), valid = self.r.wdw_eq(row, col, self._dd_us)
        return (row_us, col_us), valid

    def _pit(self, row=None, col=None):
        """pit predicate"""
        if (row is None) and (col is None):  # full raster
            return self.r == self.pit_value
        else:  # selected points
            return self.r[row, col] == self.pit_value

    # flux functions
    def _set_dd(self):
        """set drainage directions for all valid cells"""
        r0, c0 = np.where(self._mask2d)
        (self._r_ds, self._c_ds), valid = self._nb_downstream_idx(r0, c0)
        self._r, self._c = r0[valid], c0[valid]
        self._idx_ds = np.ravel_multi_index(
            (self._r_ds, self._c_ds), self.shape[-2:])
        # keep unique indices
        self._idx_ds_u, u = np.unique(self._idx_ds, return_index=True)
        self._r_ds_u, self._c_ds_u = self._r_ds[u], self._c_ds[u]

    def dd_flux(self, state):
        if not hasattr(self, '_idx_ds_u'):
            self._set_dd()
        assert state.shape == self.shape[-2:]
        state_1d = (state[self._r, self._c]).astype(
            np.float64)  # errors occur with float32
        # group at confluence of streamflows
        flux = np.zeros_like(state)
        flux_1d = np.bincount(self._idx_ds, weights=state_1d,)[self._idx_ds_u]
        assert np.isclose(state_1d.sum(), flux_1d.sum(), rtol=1e-10)
        flux[self._r_ds_u, self._c_ds_u] = flux_1d
        return flux

    def dd_accuflux(self, initial_state, maxiter=1e5):
        converged = False
        state = initial_state.copy()
        i = 0
        while (not converged) or (i <= maxiter):
            previous_state = state.copy()
            flux = self.dd_flux(previous_state)
            state = initial_state.copy() + flux
            converged = True if np.sum(state-previous_state) == 0 else False
            assert state.max() <= initial_state.sum()
            i += 1
        return state, converged

    # spatial attributes
    def index(self, x, y):
        y, x = np.atleast_1d(y), np.atleast_1d(x)
        col, row = ~self.transform * (x, y)
        return row, col

    def xy(self, row, col):
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        return self.transform * (col, row)

    # drainage network functionality
    def next_upstream(self, row=None, col=None):
        """Return indices of cells upstream from row, col.

        Note that multiple upstream cells can be found for one row, col index.

        Arguments
        ---------
        row : list, nd array
            row indices
        col : list, nd array
            column indices

        Returns
        -------
        y_out : nd array
            row indices of upstream cells
        x_out : nd array
            column indices of upstream cells
        """
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        return self._nb_upstream_idx(row, col)[0]

    def next_downstream(self, row, col, ignore_mask=False):
        """Return indices of cells downstream from row, col.

        Note that for row, col indices which are at a pit/outlet the same row, col
        indices are returned

        Arguments
        ---------
        row : list, nd array
            row indices
        col : list, nd array
            column indices

        Returns
        -------
        y_out : nd array
            row indices of downstream cells
        x_out : nd array
            column indices of downstream cells
        """
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        return self._nb_downstream_idx(row, col)[0]

    def get_pits(self, return_index=True):
        """Return the indices of outlets/pits in whole domain.

        Returns
        -------
        row : nd array
            row indices of upstream cells
        col : nd array
            column indices of upstream cells
        """
        pit_bool = self.is_pit()  # 2d array
        if return_index:
            if pit_bool.ndim == 3:
                pit_bool = pit_bool[0, :, :]
            row, col = np.where(pit_bool)
            return np.atleast_1d(row), np.atleast_1d(col)
        else:
            return pit_bool

    def is_pit(self, row=None, col=None):
        """Predicate for a pit/outlet at given indices col, row.

        Arguments
        ---------
        row : list, nd array (optional)
            row indices
        col : list, nd array (optional)
            column indices

        Returns
        -------
        pits : boolean nd array
            True where pit
        """
        if (row is not None) and (col is not None):
            row, col = np.atleast_1d(row), np.atleast_1d(col)
        return self._pit(row, col)

    def get_path_upstream(self, row, col, n=None, stop_flags=None):
        """Move upstream along drainage direction and return path.
        Stop when one of the following criteria is met:
        - no upstream cell found
        - n iterations reached (if <n> is not None)
        - at confluence
        - if it reaches a 1 (True) in <stop_flags> array

        Arguments
        ---------
        row : list, nd array
            row indices
        col : list, nd array
            column indices
        n : int (optional)
            maximum number of iteration
        stop_flags : nd array
            mask with

        Returns
        -------
        path : list
            list of (row, col) tuples
        row : int
            next upstream row index if any, empty array if none
        col : int
            next upstream column index if any, empty array if none
        """

        # init
        if (stop_flags is not None):
            if not (stop_flags.shape == self.shape[-2:]):
                msg = "The stop_flags array does not have the same shape as the drainage direction array."
                raise ValueError(msg)
            if not np.issubdtype(type(stop_flags), np.bool):
                raise ValueError(
                    "The stop_flags array is not of boolean dtype.")
        assert np.issubdtype(type(row), np.integer), np.issubdtype(
            type(col), np.integer)
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        # loop
        path, i, n_up = [], 0, 1
        while n_up >= 1:
            i += 1
            path.append((row[0], col[0]))
            row, col = self.next_upstream(row, col)
            n_up = len(row)
            # stop if flag in stop_flags or after n iterations
            if (i == n) or ((stop_flags is not None) and stop_flags[row, col]):
                path.append((row[0], col[0]))
                break
        return path, row, col

    def get_path_downstream(self, row, col, n=None, stop_flags=None):
        """Move downstream along drainage direction and return path.
        Stop when one of the following criteria is met:
        - no downstream cell found (at outlet/pit)
        - n iterations reached (if <n> is not None)
        - if it reaches a 1 (True) in <stop_flags> array

        Arguments
        ---------
        row : list, nd array
            row indices
        col : list, nd array
            column indices
        n : int (optional)
            maximum number of iteration
        stop_flags : nd array
            mask with

        Returns
        -------
        path : list
            list of (row, col) tuples
        row : int
            next downstream row index if any, empty array if none
        col : int
            next downstream column index if any, empty array if none
        """

        # init
        if (stop_flags is not None):
            if not (stop_flags.shape == self.shape[-2:]):
                msg = "The stop_flags array does not have the same shape as the drainage direction array."
                raise ValueError(msg)
            if not np.issubdtype(type(stop_flags), np.bool):
                raise ValueError(
                    "The stop_flags array is not of boolean dtype.")
        assert np.issubdtype(type(row), np.integer), np.issubdtype(
            type(col), np.integer)
        row, col = np.atleast_1d(row), np.atleast_1d(col)

        # loop
        path, i, at_outlet = [], 0, False
        while not at_outlet:
            i += 1
            path.append((row[0], col[0]))
            row, col = self.next_downstream(row, col)
            at_outlet = self.is_pit(row, col) or (len(row) == 0)
            if (i == n) or ((stop_flags is not None) and stop_flags[row, col]):
                path.append((row[0], col[0]))
                break
        return path, row, col

    def get_catchment(self, row, col, idx=None):
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        fill_value = self.r.fill_value
        catchment = np.ma.ones(self.r.shape[-2:]) * fill_value
        for (r, c, i) in zip(row, col, idx):
            while True:
                catchment[r, c] = i
                r, c = self.next_upstream(r, c)
                if r.size == 0:
                    break
        return np.ma.masked_equal(catchment, fill_value, copy=False)

    # plot
    def plot_dd(self, ax=None, **kwargs):
        import matplotlib.pyplot as plt
        if not hasattr(self, '_idx_ds_u'):
            self._set_dd()
        if ax is None:
            fig = plt.figure(figsize=self.shape[-2:])
            ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        ax.quiver(self._c, self._r, self._c_ds-self._c, self._r -
                  self._r_ds, units='xy', scale=1, **kwargs)
        ax.set_title('{} quiver'.format(self.type))
        return ax


class NextXY(DrainageDirection):
    def __init__(self, nextxy, transform,
                 search_radius=1, fill_value=-9999, zero_based=False, **kwargs):
        """initialize drainage direction object based on nextxy definition"""

        assert (nextxy.shape[0] == 2) and (
            nextxy.ndim == 3), "check nextxy dimensions"
        # convert fortran to python numbering
        if not zero_based:
            nextxy = np.where(nextxy > 0, nextxy - 1, nextxy)
        # pits definition
        pit_value = np.array([-9, -9])[:,  None, None]

        super(NextXY, self).__init__('nextxy', dd=None,
                                     raster=nextxy, transform=transform,
                                     pit_value=pit_value, fill_value=fill_value,
                                     search_radius=search_radius, **kwargs)

    # overwrite general methods in parent object
    def _nb_downstream_idx(self, row, col):
        """neighborhood downstream index"""
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        # skip outlets and outside mask
        is_valid = np.logical_not(
            self._pit(row, col), self.r.mask[0, row, col])
        col_ds, row_ds = self.r[:, row[is_valid], col[is_valid]]
        return (row_ds, col_ds), is_valid

    def _nb_upstream_idx(self, row, col):
        """neighborhood upstream index"""
        xy = np.array([col, row]).reshape(2, len(row), 1)
        (row_us, col_us), is_valid = self.r.wdw_eq(row, col, xy)
        row_us, col_us = np.atleast_1d(row_us), np.atleast_1d(col_us)
        return (row_us, col_us), is_valid

    def _pit(self, row, col):
        """pit predicate"""
        if (row is None) and (col is None):
            return self.r == self.pit_value
        else:
            return np.all(self.r[:, row, col] == self.pit_value, axis=(0, 1))


class LDD(DrainageDirection):
    def __init__(self, ldd, transform,
                 fill_value=255, **kwargs):
        """initialize drainage direction object based on pcraster ldd definition,
        for more info on the ldd data format see:
        http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/secdatbase.html#formldd
        """

        # ldd defintion
        dd = np.array([[7, 8, 9], [4, 5, 6], [1, 2, 3]])
        pit_value = 5
        super(LDD, self).__init__('ldd', dd=dd, raster=ldd, transform=transform,
                                  pit_value=pit_value, fill_value=fill_value,
                                  search_radius=1, **kwargs)


class D8(DrainageDirection):
    def __init__(self, d8, transform,
                 fill_value=255, **kwargs):
        """initialize drainage direction object based on hydroshseds d8 definition,
        for more info on the d8 data format see: https://hydrosheds.cr.usgs.gov/hydro.php
        """

        # d8 defintion
        dd = np.array([[32, 64, 128], [16, 0, 1], [8, 4, 2]])
        pit_value = 0

        super(D8, self).__init__('d8', dd=dd, raster=d8, transform=transform,
                                 pit_value=pit_value, fill_value=fill_value,
                                 search_radius=1, **kwargs)


class D8plus(DrainageDirection):
    def __init__(self, d8, transform,
                 fill_value=-9, **kwargs):
        """initialize drainage direction object based on hydroshseds d8 definition,
        extended with multiple pit values for rivermouth (0) and inland depression (-1)
        for more info on the ldd data format see: http://hydro.iis.u-tokyo.ac.jp/~yamadai/GlobalDir/
        """

        # d8 plus defintion
        dd = np.array([[32, 64, 128], [16, -1, 1], [8, 4, 2]])
        self.rivmth_value = 0  # river mouth
        pit_value = -1  # inland depression
        super(D8plus, self).__init__('d8', dd=dd, raster=d8, transform=transform,
                                     pit_value=pit_value, fill_value=fill_value,
                                     search_radius=1, **kwargs)

    # overwrite parent general method
    def _pit(self, row=None, col=None):
        """pit predicate"""
        if (row is None) and (col is None):
            return np.logical_or(self.r == self.pit_value,
                                 self.r == self.rivmth_value)
        else:
            return np.logical_or(self.r[row, col] == self.pit_value,
                                 self.r[row, col] == self.rivmth_value)

# TODO: test if faster than bincount in dd_flux function


def group(values):
    values.sort()
    dif = np.ones(values.shape, values.dtype)
    dif[1:] = np.diff(values)
    idx = np.concatenate((np.where(dif)[0], [len(values)]))
    vals = values[idx[:-1]]
    count = np.diff(idx)
    return vals, count

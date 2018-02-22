#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Dirk Eilander (contact: dirk.eilancer@vu.nl)
# Created: Nov-2017
#
# SVN auto-props
# $ID: $
# $Date: 2018-02-22 15:23:32 +0100 (Thu, 22 Feb 2018) $
# $Author: der240 $
# $Revision: 219 $
# $HeadURL: https://svn.vu.nl/repos/compound_floods/trunk/nb_ops/nb_ops/dd_ops.py $

import numpy as np
from nb_ops import MaskedWdwArray
# TODO: add io, check if sound and bbox-clip options
# TODO: test if drainge area etc function work fast compared to pcraster

class DrainageDirection(object):
    """Drainage direction numpy-based object to do fast operations on (nextx, nexty),
    LDD or D8 drainge direction networks"""

    def __init__(self, stype,
                 raster, transform, nodata=np.nan,
                 search_radius=1, **kwargs):
        self.type = stype
        raster = np.ma.masked_equal(raster, nodata, copy=False)
        self.r = MaskedWdwArray(raster, sr=search_radius)
        self.transform = transform
        self.shape = self.r.shape

    # spatial attributes
    def index(self, x, y):
        y, x = np.atleast_1d(y), np.atleast_1d(x)
        col, row = ~self.transform * (x, y)
        return row, col

    def xy(self, row, col):
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        return self.transform * (col, row)

    # drainage network functionality
    def find_upstream(self, row, col, ignore_mask=False):
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
        return self._nb_upstream_idx(row, col)

    def find_downstream(self, row, col, ignore_mask=False):
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
        return self._nb_downstream_idx(row, col)

    def find_pits(self, return_index=True):
        """Return the indices of outlets/pits in whole domain.

        Returns
        -------
        row : nd array
            row indices of upstream cells
        col : nd array
            column indices of upstream cells
        """
        pit_bool = np.logical_and(self.r == self._pit_value, ~self.r.mask)
        if return_index:
            row, col = np.where(pit_bool)[-2:]
            return np.atleast_1d(row), np.atleast_1d(col)
        else:
            return pit_bool

    def is_pit(self, row, col):
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
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        return self._pit(row, col)

    def upstream_path(self, row, col, n=None, stop_flags=None):
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
                raise ValueError("The stop_flags array is not of boolean dtype.")
        assert np.issubdtype(type(row), np.integer), np.issubdtype(type(col), np.integer)
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        # loop
        path, i, n_up = [], 0, 1
        while n_up >= 1:
            i += 1
            path.append((row[0], col[0]))
            row, col = self.find_upstream(row, col)
            n_up = len(row)
            # stop if flag in stop_flags or after n iterations
            if (i == n) or ((stop_flags is not None) and stop_flags[row,col]):
                path.append((row[0], col[0]))
                break
        return path, row, col

    def downstream_path(self, row, col, n=None, stop_flags=None):
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
                raise ValueError("The stop_flags array is not of boolean dtype.")
        assert np.issubdtype(type(row), np.integer), np.issubdtype(type(col), np.integer)
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        # loop
        path, i, at_outlet = [], 0, False
        while not at_outlet:
            i += 1
            path.append((row[0], col[0]))
            row, col = self.find_downstream(row, col)
            at_outlet = self.is_pit(row, col) or (len(row) == 0)
            if (i == n) or ((stop_flags is not None) and stop_flags[row,col]):
                path.append((row[0], col[0]))
                break
        return path, row, col

    def get_catchments(self, row, col, idx=None):
        nodata = self.r.fill_value
        catchment = np.ma.ones(self.r.shape[-2:]) * nodata
        if idx is None:
            idx = np.arange(len(row))+1
        for (r, c, i) in zip(row, col, idx):
            while True:
                catchment[r, c] = i
                r, c = self.find_upstream(r, c)
                if r.size == 0:
                    break
        return np.ma.masked_equal(catchment, nodata, copy=False)


class NextXY(DrainageDirection):
    def __init__(self, nextxy, transform,
                 search_radius=1, nodata=-9999, zero_based=False, **kwargs):
        """initialize drainage direction object based on nextxy definition"""

        assert (nextxy.shape[0] == 2) and (nextxy.ndim==3), "check nextxy dimensions"
        if not zero_based:
            nextxy = np.where(nextxy>0, nextxy - 1, nextxy)
        # pits definition
        self._pit_value = np.array([-9, -9])[:,  None, None]

        super(NextXY, self).__init__('nextxy', nextxy, transform,
                                     search_radius=search_radius, nodata=nodata,
                                     **kwargs)

    def _nb_downstream_idx(self, row, col):
        """neighborhood downstream index"""
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        # skip outlets
        is_outlet = self._pit(y_out, x_out)
        nextx, nexty = np.asarray(self.r[:, row[~is_outlet], col[~is_outlet]])
        return nexty, nextx

    def _nb_upstream_idx(self, row, col):
        """neighborhood upstream index"""
        xy = np.array([col, row]).reshape(2, len(row), 1)
        return self.r.wdw_eq(row, col, xy, return_index=True)

    def _pit(self, row, col):
        """index pit predicate"""
        return self.r[:, row, col] == self._pit_value

class LDD(DrainageDirection):
    def __init__(self, ldd, transform,
                 nodata=255, **kwargs):
        """initialize drainage direction object based on pcraster ldd definition,
        for more info on the ldd data format see:
        http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/secdatbase.html#formldd
        """

        # ldd defintion
        self.ldd0 = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]])
        self._pit_value = 5
        # create up and downstream nb patterns
        self._ldd_ds = self.ldd0.flatten()
        self._ldd_ds = np.where(self._ldd_ds==self._pit_value, 0, self._ldd_ds)
        self._ldd_us = np.flipud(self.ldd0).flatten()
        self._ldd_us = np.where(self._ldd_us==self._pit_value, 0, self._ldd_us)

        super(LDD, self).__init__('ldd', ldd, transform,
                                  search_radius=1, nodata=nodata, **kwargs)

    def _nb_downstream_idx(self, row, col):
        """neighborhood downstream index"""
        return self.r.wdw_eq(row, col, self._ldd_ds, return_index=True)

    def _nb_upstream_idx(self, row, col):
        """neighborhood upstream index"""
        return self.r.wdw_eq(row, col, self._ldd_us, return_index=True)

    def _pit(self, row=None, col=None):
        """pit predicate"""
        return self._ldd[row, col] == self._pit_value[:, None]

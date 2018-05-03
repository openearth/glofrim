#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Dirk Eilander (contact: dirk.eilancer@vu.nl)
# Created: Nov-2017
#
# objective is to get fast window functionality fully based on numpy and
# witouth having to pad the array as a good compromise between speed and memory
# we make use of masked numpy arrays to deal with missing values

import numpy as np

class MaskedWdwArray(np.ma.MaskedArray):
    """cast of numpy masked array with window functions"""
    def __new__(cls, input_array, sr=1):
        assert 2 <= input_array.ndim <= 3, "This class has only been implemented for 2d and 3d arrays"
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asanyarray(input_array).view(cls)
        # add the new attributes to the created instance
        step = np.complex(sr*2+1)
        obj._wdw_size = sr*2+1
        obj._sr = sr
        obj._drow, obj._dcol = np.mgrid[-sr:sr:step, -sr:sr:step].astype(int)
        obj._drow = obj._drow.flatten()[None, :]
        obj._dcol = obj._dcol.flatten()[None, :]
        return obj

    def wdw(self, row, col):
        """Return 2d array with all values of the windows around row, col
        flattened in the second axis and the corresponding row and col indices.

        For windows at the edges: the part ouside of the raster is masked.

        Arguments
        ---------
        row : list, nd array
            row indices
        col : list, nd array
            column indices

        Returns
        -------
        wdw : nd array
            2d array with all values of the windows aroudn row (row), col (col)
        (row_wdw, col_wdw) : tuple of nd array
            row and column indices corresponding to wdw output
        """
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        row_wdw = self._drow + row[:, None]
        col_wdw = self._dcol + col[:, None]
        edge = np.logical_or.reduce((
            (row.min()-self._sr) < 0, (row.max()+self._sr) >= self.shape[-2],
            (col.min()-self._sr) < 0, (col.max()+self._sr) >= self.shape[-1]))
        # fix indices at edges and mask out parts of window outside
        if edge:
            mask = np.logical_or(
                np.ma.masked_outside(row_wdw, 0, self.shape[-2]-1).mask,
                np.ma.masked_outside(col_wdw, 0, self.shape[-1]-1).mask)
            # make sure indices are within domain
            row_wdw0 = np.clip(row_wdw, 0, self.shape[-2]-1)
            col_wdw0 = np.clip(col_wdw, 0, self.shape[-1]-1)
            if self.ndim == 3:
                mask = np.repeat(mask[None, :, :], self.shape[0], axis=0)
                data = self[:, row_wdw0, col_wdw0]
            else:
                data = self[row_wdw0, col_wdw0]
            data = np.ma.masked_where(mask, data)
        else:
            if self.ndim == 3:
                data = np.asanyarray(self[:, row_wdw, col_wdw])
            else:
                data = np.asanyarray(self[row_wdw, col_wdw])
        return data, (row_wdw, col_wdw)

    def wdw_eq(self, row, col, value, axis=-1):
        """Check if value in window. Return a boolean array or indices that match
        the value (if return_index is True).

        Note that number of output indices can be larger or smaller than the
        input indices.

        Arguments
        ---------
        row : list, nd array
            row indices
        col : list, nd array
            column indices
        value : int, float
            value to find the in window

        Returns
        -------
        (row, col) : tuple of nd arrays
            row, col indices of upstream cells
        wdw_bool : nd array
            boolean flattend representation of windows, True where value
        """
        
        value = np.atleast_1d(value)
        wdw_data, wdw_idx = self.wdw(row, col)
        wdw_bool = np.logical_and(wdw_data.data == value, wdw_data.mask == False)
        if self.ndim == 3:
            wdw_bool = np.all(wdw_bool, axis=0)
        row = np.atleast_1d(wdw_idx[0][wdw_bool])
        col = np.atleast_1d(wdw_idx[1][wdw_bool])
        valid = wdw_bool.any(axis=-1)
        return (row, col), valid

    def wdw_max(self, *index):
        return self.wdw(*index)[0].max(axis=-1)

    def wdw_min(self, *index):
        return self.wdw(*index)[0].min(axis=-1)

    def wdw_mean(self, *index):
        return self.wdw(*index)[0].mean(axis=-1)

#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Dirk Eilander (contact: dirk.eilancer@vu.nl)
# Created: Nov-2017

from dd_ops import LDD, NextXY
import rasterio

def read_dd_rasterio(fn, ddtype='ldd', ddkwargs={}):
    with rasterio.open(fn, 'r') as src:
        if ddtype == 'ldd':
            dd_r = LDD(src.read(1), src.transform, nodoata=src.nodata, **ddkwargs)
        elif ddtype == 'nextxy':
            dd_r = NextXY(src.read(), src.transform, nodata=src.nodata, **ddkwargs)
    return dd_r 

def read_dd_pcraster(fn, transform, nodata=255):
    import pcraster as pcr
    lddmap = pcr.readmap(str(fn))
    ldd_data = pcr.pcr2numpy(lddmap, nodata)
    ldd = LDD(ldd_data, transform, nodata=nodata)
    return ldd 

def read_raster(fn):
    return None
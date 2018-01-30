#!/usr/local/bin/python
import os
from os.path import join
import rasterio
from subprocess import check_output, STDOUT, CalledProcessError

def read_mapfile(fn_map):
    with rasterio.open(fn_map, 'r') as ds:
        westin  = ds.bounds.left
        eastin  = ds.bounds.right
        northin = ds.bounds.top
        southin = ds.bounds.bottom
        res = ds.res

    if northin > southin:
        olat = 'StoN' # Data order
    else:
        olat = 'NtoS'
        northin = tmp
        northin = soutnin
        southin = tm

    return res, westin, eastin, northin, southin, olat

def get_res_str(res):
    if res[0] != res[1]:
        raise ValueError("Inconsistent spatial resolutions between lon and lat.")
    dict_res = { 0.5:'30min', 1:'1deg' }
    res_str = dict_res[res[0]]
    return res_str

def subcall(msg):
    try:
        output = check_output(msg, stderr=STDOUT, shell=True)
    except CalledProcessError as exc:
        print(exc.output)

# input
mapdir = '../PCR_Elbe/input30min/'
fn_map = join(mapdir, 'landmask_elbe_30min.map')
mkinclude = '../../cama-flood_bmi/adm/Mkinclude'

# main
# compilation
msg1 = './run-makefile.sh {}'.format(mkinclude)
print(msg1)
subcall(msg1)
# read domain info
res, westin, eastin, northin, southin, olat = read_mapfile(fn_map)
res_str =  get_res_str(res)
# generate inpmat
msg2 = './generate_inpmat {} {} {} {} {} {}'.format(res[0], westin, eastin, northin, southin, olat)
print(msg2)
subcall(msg2)

# output
fn_diminfo   = 'diminfo_{}.txt'.format(res_str)
fn_inpmatbin = 'inpmat-{}.bin'.format(res_str)
fn_inpmattxt = 'inpmat-{}.txt'.format(res_str)

os.rename('diminfo_tmp.txt', fn_diminfo  )
os.rename('inpmat-tmp.bin',  fn_inpmatbin)
os.rename('inpmat-tmp.txt',  fn_inpmattxt)

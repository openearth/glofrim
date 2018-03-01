import numpy as np
import os
import rasterio
import click
import glob

@click.command()
@click.argument('in_fn')
                # , required=True,
                # help='filename of binary map control (*.ctl) file or directory with control files')
@click.option('-o', '--fn_out', default=None,
                help='filename of output file name, defaults to same basename as input file')
@click.option('--epsg', default=4326, help='crs epsg code, 4326 by default.')
@click.option('--force_dtype', help='overwrite ctl file based dtype (numpy dtype format)')
@click.option('--force_endian', type=click.Choice(['little', 'big']),
                help='overwrite little / big endian setting from ctl file')
@click.option('--force_overwrite', help='force overwrite of existing output geotiff file', is_flag=True)
@click.option('--cell_center/--cell_corner', default=True)
def bin_2_tif(in_fn, fn_out=None, epsg=4326, force_overwrite=False,
              force_dtype=None, force_endian=None, cell_center=True):
    """read binary map based on control file and convert to geotiff
    """
    fn = os.path.abspath(in_fn)
    if fn_out is not None:
        fn_out = os.path.abspath(fn_out)
    fns = []
    if os.path.isdir(fn):
        fns = [fn for fn in glob.glob(os.path.join(fn, '*.ctl'))
                    if not os.path.basename(fn).startswith('inpmat')]
        if (fn_out is not None) and (len(fns) > 1):
            click.echo('fn_out is ignored as multiple output files are found')
            fn_out = None
    elif os.path.isfile(fn):
        fns = [fn]
    else:
        fns = glob.glob(fn)
    if len(fns) == 0:
        raise IOError("No ctl files found matching {:s}".format(fn))

    for fn in fns:
        if fn_out is None:
            fno = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace('.ctl', '.tif'))
        else:
            fno = fn_out

        if force_overwrite and os.path.isfile(fno):
            os.unlink(fno)
        elif os.path.isfile(fno):
            click.echo("File {} already exists. Use --force_overwrite to overwrite.".format(os.path.basename(fno)))
            continue

        click.echo("processing {:s}".format(os.path.basename(fn)))
        data, prof = read_bin(fn, dtype=force_dtype, force_endian=force_endian,
                                cell_center=cell_center)
        # write to geotiff
        crs = {'init': 'epsg:{:d}'.format(int(epsg))}
        prof.update(driver='GTiff', compress='lzw', dtype=str(data.dtype), crs=crs)
        with rasterio.open(fno, 'w', **prof) as dst:
            for i in xrange(prof['count']):
                dst.write(data[i, :, :], i+1)

def read_bin(fn_ctl, dtype=None, force_endian=None, cell_center=True):
    # read data props
    ctl = ctl_parser(fn_ctl)
    # convert to rasterio profile taxonomy
    prof = ctl_2_rasterio_profile(ctl, dtype=dtype, cell_center=cell_center)
    # read data
    # import pdb; pdb.set_trace()
    nlayer, nrow, ncol = prof['count'], prof['height'], prof['width']
    if force_endian is not None:
        if force_endian not in ['big', 'little']:
            raise ValueError('force_endian must be "big" or "little"')
        force_endian = '>' if force_endian=='big' else '<'
        prof['dtype'] = force_endian + prof['dtype'][1:]
    try:
        data = np.fromfile(ctl['fn'], prof['dtype']).reshape(nlayer, nrow, ncol)
    except ValueError:
        if 'f4' in prof['dtype']:
            prof['dtype'] = prof['dtype'].replace('f4', 'f8')
        data = np.fromfile(ctl['fn'], prof['dtype']).reshape(nlayer, nrow, ncol)
    if not ctl['options']['yrev']:
        data = data[:, ::-1, :]
    return data, prof

# read files via ctl
def ctl_parser(fn):
    """ parse grads ctl file, return dict with all values

    for more information about the ctl files see:
    http://cola.gmu.edu/grads/gadoc/descriptorfile.html#XDEF"""
    ctl_map = {'dset': {'fmt': [str]},
        'undef': {'fmt': [int]},
        'title': {'fmt': [str]},
        'options': {'fmt': [str, str], 'names': ['yrev', 'endian']},
        'xdef': {'fmt': [int, str, float, float],
                 'names': ['n', 'mapping', 'start', 'increment']},
        'ydef': {'fmt': [int, str, float, float],
                 'names': ['n', 'mapping', 'start', 'increment']},
        'zdef': {'fmt': [int, str, float, float],
                 'names': ['n', 'mapping', 'start', 'increment']},
        'tdef': {'fmt': [int, str, str, str],
                 'names': ['n', 'mapping', 'start', 'increment']},
        'vars': {'fmt': [int]},
    }
    # mapping from grads to numpy dtype
    vars_unity_dtype = {'-1,40,2,-1': 'i2',
                        '-1,40,4': 'i4',
                        '99': 'f4'}
    var_map, var_fmt = ['nlayers', 'units'], [int, str]

    with open(fn, 'r') as f:
        lines = f.readlines()

    # read definitions
    ctl = {}
    for name in ctl_map:
        fmt = ctl_map[name]['fmt']
        if len(fmt) == 1:
            ctl[name] = get_ctlvalue(name, lines, fmt)[0]
        else:
            names = ctl_map[name]['names']
            ctl[name] = dict(zip(names, get_ctlvalue(name, lines, fmt)[0]))
    # transform into boolean
    ctl['options']['yrev'] = ctl['options']['yrev'] == 'yrev'
    # strip fn
    if ctl['dset'].startswith('^'):
        ctl['dset'] = ctl['dset'][1:]
    # check if known mapping types. for now just linear
    for name in ['xdef', 'ydef', 'zdef', 'tdef']:
        assert ctl[name]['mapping'] == 'linear'
    # for now ignore time and z dimensions
    assert ctl['tdef']['n'] == 1

    # read n variables metadata
    ctl['endian'] = '<' if ctl['options']['endian'].startswith('little') else '>'
    ivar = [i for i, line in enumerate(lines) if line.startswith('var')][0] + 1
    nvar = ctl['vars']
    ctl['vars'] = {}
    var_names = []
    for i in xrange(nvar):
        var_name = lines[ivar + i].split()[0].strip()
        var_names.append(var_name)
        values, comments = get_ctlvalue(var_name, lines, var_fmt)
        ctl['vars'][var_name] = dict(zip(var_map, values))
        ctl['vars'][var_name]['ilayer'] = i
        ctl['vars'][var_name]['description'] = comments

        # check for known datatypes others should be implemented, see link above for more datatype
        if ctl['vars'][var_name]['units'] in vars_unity_dtype.keys():
            # map grads dtype to numpy dtype
            ctl['vars'][var_name]['units'] = vars_unity_dtype[ctl['vars'][var_name]['units']]
        else:
            raise ValueError("unknown dtype")

    ctl['var_names'] = var_names
    ctl['fn'] = os.path.join(os.path.dirname(fn), ctl['dset'])
    return ctl

def get_ctlvalue(name, lines, fmt, ignore='**'):
    """function to read values from a *.ctl file"""
    values = []
    n = len(fmt)
    for line in lines:
        if line.startswith(name + ' '):
            linesplit = line.split(ignore)
            settings = linesplit[0]
            if len(linesplit) > 1:
                comments = ignore.join(linesplit[1:]).strip()
            else:
                comments = ''
            settings = settings.replace(name, '')
            if n > 1:
                v = settings.split()
                values = [f(v[i].strip()) for i,f in enumerate(fmt)]
            else:
                values = fmt[0](settings.strip())
            break
    return values, comments


def ctl_2_rasterio_profile(ctl, dtype=None, cell_center=True):
    from rasterio.transform import from_origin
    prof = {}
    if dtype is None:
        prof['dtype'] = ctl['endian'] + ctl['vars'][ctl['var_names'][-1]]['units'] #assume same dtype for all variables
    else:
        prof['dtype'] = dtype
    prof['width'], prof['height'] = ctl['xdef']['n'], ctl['ydef']['n']
    prof['count'] = np.sum([ctl['vars'][name]['nlayers'] for name in ctl['var_names']])
    prof['nodata'] = ctl['undef']
    # calculate transform
    resx = float(np.abs(ctl['xdef']['increment']))
    resy = float(np.abs(ctl['ydef']['increment']))
    xmin = ctl['xdef']['start']
    if ctl['options']['yrev']:
        ymin = ctl['ydef']['start']
        ymax = ymin + ctl['ydef']['n']*ctl['ydef']['increment']
    else:
        ymax = ctl['ydef']['start']
    if cell_center:
        prof['transform'] = from_origin(xmin-resx/2., ymax-resy/2., resx, resy)
    else:
        prof['transform'] = from_origin(xmin, ymax, resx, resy)
    return prof

if __name__ == "__main__":
    bin_2_tif()

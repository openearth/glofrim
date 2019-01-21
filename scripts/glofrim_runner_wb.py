from glofrim import Glofrim 

import click
from os.path import isdir, dirname, abspath
from os import makedirs
import os
from dateutil import parser
import pandas as pd
import numpy as np

# utils
def parse_datetime(param, value):
    try:
        date = parser.parse(value)
    except:
        raise click.BadParameter("Couldn't understand date for the '{}' argument.".format(param))
    return date

def parse_dir(param, path):
    try:
        path = abspath(path)
        if not isdir(path):
            os.makedirs(path)
    except:
        raise click.BadParameter("Couldn't understand or create folder directory for the '{}' argument.".format(param))
    return path

@click.group()
def cli():
    pass

@click.command()
@click.argument('ini',)
@click.option('--env', default='', help='path to glofrim env file with engine paths')
@click.option('-o', '--out-dir', default='', help='directory to save model outputs', type=click.Path())
@click.option('-s', '--start-date', default='', help='set start time for all models')
@click.option('-e', '--end-date', default='', help='set end time for all models')
def run(ini, env='', out_dir='', end_date='', start_date=''):
    """
    Runs a coupled model using the GLOFRIM API.

    The model run times are based on the individual model configuration files unless set with the 'start-date' and 'end-date' arguments.
    The coupled model starts at the latest start date of all individual models. 
    The individual models will first be updated until the combined start date before any exchage between the models takes place.
    Any forcing that is not provided from a coupled model via BMI should be set in the individual configuration files.     

    The ini argument contains information about the engine (if the model is not a python package), and path to the configuration file the individual models, 
    as well as the timestep and variables of exchange between the models. An example is given of the glofrim.ini file is given in the glofrim root dir.

    EXAMPLE

    python glofrim_runner.py run /path/to/glofrim.ini --env /path/to/glofrim.env -s 200-01-01 -e 2001-01-01
    """
    # validate and parse arguments
    if end_date:
        end_date = parse_datetime('end-date', end_date)
    if start_date:
        start_date = parse_datetime('start-date', start_date)
    if start_date and end_date:
        if end_date <= start_date:
            raise click.BadParameter("'end-time' should be larger than start-time")
    if out_dir:
        out_dir = parse_dir('out-dir', out_dir)
    # initialize combined bmi
    cbmi = Glofrim()
    # initialize configuration
    ini = abspath(ini)
    print(ini)
    cbmi.initialize_config(config_fn=ini, env_fn=env)
    # get model start and end times 
    if start_date:
        cbmi.set_start_time(start_date)
    if end_date:
        cbmi.set_end_time(end_date)
    # create and set optional outdir 
    if isdir(out_dir):
        cbmi.set_out_dir(out_dir)
    # initialize model
    cbmi.initialize_model()
    end = cbmi.get_end_time()
    start = cbmi.get_start_time()

    # get outlets for CMF
    if 'CMF' in cbmi.bmimodels:
        # get outlets
        row, col = cbmi.bmimodels['CMF'].grid._dd.get_pits()
        pit_inds = cbmi.bmimodels['CMF'].grid.ravel_multi_index(row, col)
        # get outlet coordinates from 1d network
        coords = np.asarray(cbmi.bmimodels['CMF'].grid._1d.nodes)
        inds = cbmi.bmimodels['CMF'].grid._1d.inds
        idx = np.array([[i for i,ind in enumerate(inds) if ind == pit_ind] for pit_ind in pit_inds]).squeeze()
        pit_coords = coords[idx]
        # create dataframe
        df_pits = pd.DataFrame({'lon': pit_coords[:,0], 'lat': pit_coords[:,1], 'flat_ind': pit_inds}).set_index('flat_ind')
        if 'LFP' in cbmi.bmimodels:
            df_pits['in_lfp_domain'] = cbmi.bmimodels['LFP'].grid._inside_bounds(pit_coords[:,0], pit_coords[:,1])
        print(ini.replace('.ini', '_CMFrivmth_loc.csv'))
        df_pits.to_csv(ini.replace('.ini', '_CMFrivmth_loc.csv'))
        df_outflw = pd.DataFrame(index=pd.date_range(start, end, freq=cbmi._dt), columns=df_pits.index)
        df_outflw.index.name = 'date'
    
    # run all models until combined start time
    for mod in cbmi.bmimodels:
        if cbmi.bmimodels[mod].get_current_time() < start:
            cbmi.bmimodels[mod].update_until(start)

    # run combined models until combined end time
    # cbmi.update_until(end)
    # update per timestep to read CMF discharge outputs
    while cbmi.get_current_time() < end:
        t = cbmi.get_current_time()
        cbmi.update()

        # read discharge at outputs
        if 'CMF' in cbmi.bmimodels:
            df_outflw.loc[t, df_pits.index] = cbmi.get_value_at_indices('CMF.outflw', df_pits.index.values)

    # run indiviudal models until end
    for mod in cbmi.bmimodels:
        if cbmi.bmimodels[mod].get_current_time() < end:
            cbmi.bmimodels[mod].update_until(end)

    # finalize
    cbmi.finalize()

    if 'CMF' in cbmi.bmimodels:
        df_outflw.to_csv(ini.replace('.ini', '_CMFrivmth_outflw.csv'))

cli.add_command(run)


if __name__ == '__main__':
    cli()
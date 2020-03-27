from glofrim import Glofrim

import click
from os.path import isdir, dirname, abspath
from os import makedirs
import os
from dateutil import parser

# utils


def parse_datetime(param, value):
    try:
        date = parser.parse(value)
    except:
        raise click.BadParameter(
            "Couldn't understand date for the '{}' argument.".format(param))
    return date


def parse_dir(param, path):
    try:
        path = abspath(path)
        if not isdir(path):
            os.makedirs(path)
    except:
        raise click.BadParameter(
            "Couldn't understand or create folder directory for the '{}' argument.".format(param))
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
            raise click.BadParameter(
                "'end-time' should be larger than start-time")
    if out_dir:
        out_dir = parse_dir('out-dir', out_dir)
    # initialize combined bmi
    cbmi = Glofrim()
    # initialize configuration
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
    # run all models until combined start time
    for mod in cbmi.bmimodels:
        if cbmi.bmimodels[mod].get_current_time() < start:
            cbmi.bmimodels[mod].update_until(start)
    # run combined models until combined end time
    cbmi.update_until(end)
    # run indiviudal models until end
    for mod in cbmi.bmimodels:
        if cbmi.bmimodels[mod].get_current_time() < end:
            cbmi.bmimodels[mod].update_until(end)
    # finalize
    cbmi.finalize()


@click.command()
@click.argument('ini')
@click.argument('mod')
@click.option('--env', default='', help='path to glofrim env file with engine paths')
@click.option('-o', '--out-dir', default='', help='directory to save model outputs', type=click.Path())
@click.option('-s', '--start-date', default='', help='set start time for all models')
@click.option('-e', '--end-date', default='', help='set end time for all models')
def run_single(ini, mod, env='', out_dir='', end_date='', start_date=''):
    """
    Runs a single model using the GLOFRIM API.

    The model run times are based on the model configuration file unless set with the 'start-date' and 'end-date' arguments.
    Forcing should be set in the model configuration file.     

    The ini argument contains information about the engine (if the model is not a python package), and path to the configuration file the individual models, 
    as well as the timestep and variables of exchange between the models. An example is given of the glofrim.ini file is given in the glofrim root dir.

    EXAMPLE

    python glofrim_runner.py run_singe /path/to/glofrim.ini ModelName
    """
    # validate and parse arguments
    if end_date:
        end_date = parse_datetime('end-date', end_date)
    if start_date:
        start_date = parse_datetime('start-date', start_date)
    if start_date and end_date:
        if end_date <= start_date:
            raise click.BadParameter(
                "'end-time' should be larger than start-time")
    if out_dir:
        out_dir = parse_dir('out-dir', out_dir)
    # initialize bmi
    cbmi = Glofrim()
    # initialize configuration
    cbmi.initialize_config(config_fn=ini, env_fn=env)
    if mod not in cbmi.bmimodels.keys():
        raise click.BadParameter(
            "Model with name {} not found in ini file".format(mod))
    bmi = cbmi.bmimodels[mod]
    # set model start and end time
    if start_date:
        bmi.set_start_time(start_date)
    if end_date:
        bmi.set_end_time(end_date)
    # create and set optional outdir
    if isdir(out_dir):
        bmi.set_out_dir(out_dir)
    # initialize model
    cbmi.logger.info('initializing {} model'.format(mod))
    bmi.initialize_model()
    start = bmi.get_start_time()
    end = bmi.get_end_time()
    # run combined model until end time
    cbmi.logger.info(
        'running model from {} to {}'.format(str(start), str(end)))
    bmi.update_until(end)
    # finalize
    cbmi.logger.info('finalizing model')
    cbmi.finalize()

# @click.command()
# @click.argument('ini')
# @click.argument('out_fn')
# @click.option('-e', '--env_fn', default='', help='path to glofrim env file with engine pahts')
# def spatial(ini, out_fn, env_fn=''):
#     """
#     Returns a json file with the coupled indices of two models
#     """
#     # initialize combined bmi
#     cbmi = Glofrim()
#     # initialize configuration
#     cbmi.initialize_config(config_fn=ini, env_fn=env_fn)
#     exchanges = cbmi.set_exchanges()
#     exchanges = [ex[1]['coupling'] for ex in exchanges if ex == 'exchange']
#     assert len(exchanges) == 1
#     for sc in exchanges:
#         sc.to_file(out_fn)


cli.add_command(run)
cli.add_command(run_single)
# cli.add_command(spatial)


if __name__ == '__main__':
    cli()

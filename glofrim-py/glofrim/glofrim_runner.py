import click
from os.path import isdir, dirname, abspath
from os import makedirs
from glofrim_bmi import Glofrim

@click.group()
def cli():
    pass

@click.command()
@click.argument('ini')
@click.option('-o', '--out_dir', default='', help='directory to save model outputs')
@click.option('-e', '--env_fn', default='', help='path to glofrim env file with engine paths')
def run(ini, env_fn='', out_dir=''):
    """
    Runs a coupled hydrlogical - hydrodynamic model
    """
    # initialize combined bmi
    cbmi = Glofrim()
    # initialize configuration
    cbmi.initialize_config(config_fn=ini, env_fn=env_fn)
    # create and set optional outdir 
    if (out_dir is not None) and (out_dir != ''):
        out_dir = abspath(out_dir)
        if not isdir(out_dir):
            os.makedirs(out_dir)
        cbmi.set_out_dir(out_dir)
    # initialize model
    cbmi.initialize_model()
    # run until end
    cbmi.update_until(cbmi.get_end_time())
    # finalize
    cbmi.finalize()

@click.command()
@click.argument('ini')
@click.argument('out_fn')
@click.option('-e', '--env_fn', default='', help='path to glofrim env file with engine pahts')
def spatial(ini, out_fn, env_fn=''):
    """
    Returns a json file with the coupled indices of two models
    """
    # initialize combined bmi
    cbmi = Glofrim()
    # initialize configuration
    cbmi.initialize_config(config_fn=ini, env_fn=env_fn)
    exchanges = cbmi.set_exchanges()
    exchanges = [ex[1]['coupling'] for ex in exchanges if ex == 'exchange']
    assert len(exchanges) == 1
    for sc in exchanges:
        sc.to_file(out_fn)

cli.add_command(run)
cli.add_command(spatial)

if __name__ == '__main__':
    cli()

__author__ = "DirkEilander"

import functools
import glofrim
import glofrim.glofrim_lib as glib
import os
from os.path import join, dirname, realpath, isfile, relpath, isdir
import sys
import glob
from datetime import datetime, timedelta
from termcolor import colored
import shutil
import numpy as np

ddir = dirname(realpath(__file__))

# wrapper to check if function raises error


def tryexcept(msg, errors=(Exception, ), fail=False):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            prefix = '{} - '.format(self.prefix)
            try:
                out = func(self, *args, **kwargs)
                print(prefix + msg + colored('PASS', 'green'))
                return out
            except errors, e:
                print(prefix + msg + colored('FAILED', 'red'))
                self.errors += 1
                if fail:
                    raise errors(e)
                else:
                    print(colored('ERROR: ', 'red') + repr(e))
                    return None
        return wrapper
    return decorator


class integration_test(object):

    def __init__(self, mod):
        self.mod = mod
        self.prefix = mod
        self.errors = 0
        self.bmi = None
        self.glofrim_config_fn = config_fn = join(
            ddir, r'integration_test.ini')
        # do test
        if self.mod == 'GLOFRIM':
            self.test_coupled()
        else:
            self.test_single()

    @tryexcept('set outdir: ', fail=True)
    def create_output_dir(self):
        self.outputDir = join(ddir, 'model_test_data', 'output', self.mod)
        if not isdir(self.outputDir):
            os.makedirs(self.outputDir)

    def cleanup(self):
        if isdir(self.outputDir):
            shutil.rmtree(dirname(self.outputDir))
        if self.mod != 'WFL':
            if self.bmi.initialized:
                os.unlink(self.bmi._config_fn)
        if self.mod == 'CMF':
            if hasattr(self.bmi, '_mapdir'):
                for fn in glob.glob(join(self.bmi._mapdir, '*tmp*')):
                    os.unlink(fn)
                ext = ("*.txt", "*.log", "*.nc",)
                for e in ext:
                    for fn in glob.glob(join(dirname(self.bmi._config_fn), e)):
                        os.unlink(fn)

    @tryexcept('get BMI attrs: ', fail=True)
    def get_bmi_attrs(self):
        if self.mod == 'GLOFRIM':  # coupled run
            self.config_fn = self.glofrim_config_fn

        else:  # single run
            if not isfile(self.glofrim_config_fn):
                raise IOError('glofrim inifile not found: {}'.format(
                    self.glofrim_config_fn))
            _config = glib.configread(self.glofrim_config_fn)

            if _config.has_option('models', self.mod):
                config_fn = _config.get('models', self.mod)
                self.config_fn = str(
                    join(ddir, 'model_test_data', _config.get('models', self.mod)))
            else:
                msg = 'config filename not found {} model'
                raise ValueError(msg.format(self.mod))

            # check if bmi component requires engine
            self.requires_engine = 'engine' in getattr(
                glofrim, self.mod).__init__.__code__.co_varnames
            if self.requires_engine:
                if _config.has_option('engines', self.mod):
                    self.engine = glib.getabspath(
                        _config.get('engines', self.mod), ddir)
                else:
                    msg = 'engine path not found for {} model'
                    raise ValueError(msg.format(self.mod))

    @tryexcept('initialize BMI: ', fail=True)
    def test_init_standalone_bmi(self, _loglevel='INFO'):
        bmi = getattr(glofrim, self.mod)
        bmi_kwargs = {'loglevel': _loglevel}
        if self.requires_engine:
            bmi_kwargs.update(engine=self.engine)
        return bmi(**bmi_kwargs)

    @tryexcept('initialize BMI: ', fail=True)
    def test_init_coupled_bmi(self, _loglevel='INFO'):
        bmi = getattr(glofrim, self.mod)
        bmi_kwargs = {'loglevel': _loglevel, config_fn = self.config_fn}
        return bmi(**bmi_kwargs)

    @tryexcept('initialize config: ', fail=True)
    def test_initialize_config(self):
        print(str(self.config_fn))
        self.bmi.initialize_config(str(self.config_fn))

    @tryexcept('set model specific settings: ', fail=True)
    def test_edit_config(self):
        if self.mod == 'PCR':
            inputDir = join(ddir, 'model_test_data', 'PCR_Elbe', 'input30min')
            self.bmi.set_attribute_value('globalOptions:inputDir', inputDir)
            self.bmi.set_attribute_value(
                'globalOptions:outputDir', self.outputDir)
        elif self.mod == 'CMF':
            self.bmi.set_attribute_value('OUTPUT:COUTDIR', relpath(
                self.outputDir, dirname(self.config_fn)))

    @tryexcept('set CMF inpmat: ')
    def test_set_inpmat(self):
        self.bmi.set_inpmat((7.0, 48.0, 17.0, 55.0), 0.5)
        self.bmi.set_inpmat_attrs()

    @tryexcept('set simulation time: ')
    def test_set_sim_times(self, start_time, end_time):
        self.bmi.set_start_time(start_time)
        self.bmi.set_end_time(end_time)

    @tryexcept('initialize model: ', fail=True)
    def test_initialize_model(self):
        self.bmi.initialize_model()

    @tryexcept('check simulation time: ')
    def test_assert_sim_times(self, start_time, end_time):
        if not (end_time == self.bmi.get_end_time()):
            raise AssertionError('end time not set correctly')
        if not (start_time == self.bmi.get_start_time()):
            raise AssertionError('start time not set correctly')

    @tryexcept('check datatype meta functions: ')
    def test_types(self):
        _check = [
            isinstance(self.bmi.get_component_name(), str),
            isinstance(self.bmi.get_input_var_names(), list),
            isinstance(self.bmi.get_output_var_names(), list)
        ]
        assert np.all(_check), "invalid output types from meta functions"

    @tryexcept('get grid: ')
    def test_get_grid(self):
        self.bmi.get_grid()

    @tryexcept('grid index: ')
    def test_grid_index(self):
        # import pdb; pdb.set_trace()
        if self.bmi.grid.mask is not None:
            idx0 = np.where(self.bmi.grid.mask.flat)[0][0]
        else:
            idx0 = 0
        x, y = self.bmi.grid.xy(idx0)
        idx = self.bmi.grid.index(x, y)[0]
        if self.bmi.grid.type == 1:
            # NOTE: for ucgrid the x, y coordinates can return a unit catchment located at another index
            assert idx == idx0, "Model index and xy functions do not align at set index"

    @tryexcept('model update: ')
    def test_update(self):
        try:
            # update for one second
            self.bmi.update((self.bmi._dt/2).total_seconds())
            if self.mod != 'LFP':
                raise AssertionError(
                    'Updates smaller than timestep should yield an assertion error')  # except for LFP
        except ValueError:
            pass  # the above should raise a value error
        # update for one days
        self.bmi.update(1*86400)  # this should work

    @tryexcept('model update until: ')
    def test_update_until(self):
        next_t = self.bmi.get_current_time() + (self.bmi._dt * 2)
        self.bmi.update_until(next_t)

    @tryexcept('get and set values: ')
    def test_get_set_var(self):
        var_name = self.bmi.get_output_var_names()[0]
        var_shape = self.bmi.get_var_shape(var_name)
        self.bmi.set_value(var_name, np.ones(var_shape))
        if not self.bmi.get_value(var_name).sum() == np.ones(var_shape).sum():
            raise AssertionError(
                'Get and set values does not yield the same sum')

    @tryexcept('get and set values at indices: ')
    def test_get_set_var_at_indices(self):
        var_name = self.bmi.get_output_var_names()[0]
        if self.bmi.grid._1d is not None:
            inds = self.bmi.grid._1d.inds
        else:
            inds = [0]
        self.bmi.set_value_at_indices(var_name, inds, np.ones_like(inds))
        if not (self.bmi.get_value_at_indices(var_name, inds).sum() == np.ones_like(inds).sum()):
            raise AssertionError(
                'Get and set values does not yield the same sum')

    @tryexcept('finalize model: ')
    def test_finalize(self):
        self.bmi.finalize()

    @tryexcept('stand alone BMI test: ')
    def test_single(self):
        try:
            # create output dir & read glofrim ini
            self.create_output_dir()
            self.get_bmi_attrs()

            # initialize bmi
            self.bmi = self.test_init_standalone_bmi()

            # initialize config
            self.test_initialize_config()

            # set model specific settings
            self.test_edit_config()

            # test CMF inpmap
            if self.mod == 'CMF':
                self.test_set_inpmat()

            # test settingn sim times
            start_time = datetime(2000, 2, 1)
            end_time = datetime(2000, 4, 1)
            self.test_set_sim_times(start_time, end_time)
            self.test_assert_sim_times(start_time, end_time)

            # initialize model
            self.test_initialize_model()

            # check if start & end time still the same after initialization
            self.test_assert_sim_times(start_time, end_time)

            # # check some attribute types
            self.test_types()

            # get grid
            self.test_get_grid()
            if self.bmi.grid.type != 2:
                self.test_grid_index()

            # test updates
            self.test_update()
            self.test_update_until()

            # check get_var / set_var (at indices)
            self.test_get_set_var()
            self.test_get_set_var_at_indices()

            # finalize
            self.test_finalize()

        # cleanup
        finally:
            if self.errors != 0:
                print(self.mod + " - test result: " +
                      colored("{:d} errors".format(self.errors), 'red'))
                raise AssertionError('Errors found in test')
            else:
                print(self.mod + " - test result: " +
                      colored("zero errors", 'green'))
            self.cleanup()

    @tryexcept('coupled BMI test: ')
    def test_single(self):
        try:
            # create output dir & read glofrim ini
            self.create_output_dir()
            self.get_bmi_attrs()

            # initialize bmi
            self.bmi = self.test_init_standalone_bmi()

            # initialize config
            self.test_init_coupled_bmi()

            # test settingn sim times
            start_time = datetime(2000, 2, 1)
            end_time = datetime(2000, 4, 1)
            self.test_set_sim_times(start_time, end_time)
            self.test_assert_sim_times(start_time, end_time)

            # initialize model
            self.test_initialize_model()

            # check if start & end time still the same after initialization
            self.test_assert_sim_times(start_time, end_time)

            # test updates
            self.test_update()
            self.test_update_until(end_time)

            # check get_var / set_var (at indices)
            self.test_get_set_var()
            self.test_get_set_var_at_indices()

            # finalize
            self.test_finalize()

        # cleanup
        finally:
            if self.errors != 0:
                print(self.mod + " - test result: " +
                      colored("{:d} errors".format(self.errors), 'red'))
                raise AssertionError('Errors found in test')
            else:
                print(self.mod + " - test result: " +
                      colored("zero errors", 'green'))
            self.cleanup()


if __name__ == "__main__":
    integration_test(sys.argv[1])

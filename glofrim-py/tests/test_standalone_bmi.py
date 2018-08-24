
__author__ = "DirkEilander"

import unittest
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
# env_fn = glob.glob(join(ddir, r'environment.env')) + \
#           glob.glob(join(ddir, r'../../environment.env'))
# env_fn = env_fn[0]        

class MyTest(unittest.TestCase):
    mod = 'PCR'
    env_fn = join(ddir, r'environment.env')
    bmi = None
    outputDir = None
    
    def get_config_fn(self):
        config_fn = join(ddir, r'glofrim_standalone.ini')
        config = glib.configread(config_fn)
        config_fn = str(join(ddir, 'model_test_data', config.get('models', self.mod)))
        if not isfile(config_fn):
            raise AssertionError('{} not found'.format(config_fn))
        return config_fn

    def init_standalone_bmi(self, _loglevel='INFO'):
        self.bmi = getattr(glofrim, self.mod)
        bmi_kwargs = {'loglevel': _loglevel}
        # check if bmi component requires engine
        if 'engine' in self.bmi.__init__.__code__.co_varnames:
            config = glib.configread(self.env_fn)
            if not config.has_option('engines', self.mod):
                msg = 'GLOFRIM ini or environment file misses a "engines" section with {} option'
                raise ValueError(msg.format(self.mod))
            engine_path = glib.getabspath(config.get('engines', self.mod), dirname(self.env_fn))
            bmi_kwargs.update(engine = engine_path)
        return self.bmi(**bmi_kwargs)

    def test_init_step1(self):
        try:
            self.bmi = self.init_standalone_bmi()
            print('{} - initialize BMI: '.format(self.mod) + colored('PASS', 'green'))
            
            config_fn = self.get_config_fn()
            print('{} - localize config file: '.format(self.mod))
            
            self.bmi.initialize_config(str(config_fn))
            print('{} - initialize config: '.format(self.mod) + colored('PASS', 'green'))
            
            if self.mod == 'PCR':
                inputDir = join(ddir, 'model_test_data', 'PCR_Elbe', 'input30min')
                self.outputDir = join(ddir, 'model_test_data', 'PCR_Elbe', 'output')
                self.bmi.set_attribute_value('globalOptions:inputDir', inputDir)
                self.bmi.set_attribute_value('globalOptions:outputDir', self.outputDir)
            elif self.mod == 'CMF':
                self.outputDir = join(ddir, 'model_test_data', 'CMF_Elbe', 'output')
                if not isdir(self.outputDir): os.makedirs(self.outputDir)
                self.bmi.set_attribute_value('OUTPUT:COUTDIR', relpath(self.outputDir, dirname(config_fn)))
            print('{} - change config: '.format(self.mod) + colored('PASS', 'green'))
            
            if self.mod == 'CMF':
                self.bmi.set_inpmat((7.0, 48.0, 17.0, 55.0), 0.5)
                self.bmi.set_inpmat_attrs()
                print('{} - set inpmat: '.format(self.mod) + colored('PASS', 'green'))

            start_time = datetime(2000, 2, 1)
            end_time = datetime(2000, 4, 1)
            self.bmi.set_start_time(start_time)
            self.bmi.set_end_time(end_time)
            self.assertEqual(end_time, self.bmi.get_end_time())
            self.assertEqual(start_time, self.bmi.get_start_time())
            print('{} -set/get model simulation times: '.format(self.mod) + colored('PASS', 'green'))
            
            self.bmi.initialize_model()
            print('{} - initialize model: '.format(self.mod) + colored('PASS', 'green'))
            
            # check if start & end time still the same after initialization
            self.assertEqual(end_time, self.bmi.get_end_time())
            self.assertEqual(start_time, self.bmi.get_start_time())
            print('{} - internal model simulation time: '.format(self.mod) + colored('PASS', 'green'))

            # check some attribute types
            self.assertIsInstance(self.bmi.get_component_name(), str)
            self.assertIsInstance(self.bmi.get_input_var_names(), list)
            self.assertIsInstance(self.bmi.get_output_var_names(), list)
            print('{} - attribute types: '.format(self.mod) + colored('PASS', 'green'))

            # check grid
            # if self.mod == 'CMF':
            #     glib.subcall('python ../scripts/cama_maps_io.py --force_overwrite "{}"'.format(join(self.bmi._mapdir, 'hires')))
            self.bmi.get_grid()
            print('{} - get grid info: '.format(self.mod) + colored('PASS', 'green'))
            x, y = 10.1,53.4
            self.bmi.get_grid()
            if self.bmi.grid.type != 2:
                res = self.bmi.grid.res
                idx = self.bmi.grid.index(x, y)
                _x, _y = self.bmi.grid.xy(idx)
            self.assertTrue(x-_x < res)
            self.assertTrue(y-_y < res)
            print('{} - grid index: '.format(self.mod) + colored('PASS', 'green'))

            with self.assertRaises(ValueError):
                # update for one second
                self.bmi.update(1) # this should raise a value error
            # update for one days
            self.bmi.update(1*86400) # this should work
            print('{} - model update: '.format(self.mod) + colored('PASS', 'green'))
            
            next_t = self.bmi.get_current_time() + timedelta(days=2)
            self.bmi.update_until(next_t)
            print('{} - model update_until: '.format(self.mod) + colored('PASS', 'green'))

            # check if internal model time is correct
            self.assertEquals(self.bmi._t, self.bmi.get_current_time())
            print('{} - model internal time: '.format(self.mod) + colored('PASS', 'green'))

            # check get_var / set_var
            var_name = self.bmi.get_output_var_names()[0]
            var_shape = self.bmi.get_var_shape(var_name)
            self.bmi.set_value(var_name, np.ones(var_shape)) 
            self.assertEqual(self.bmi.get_value(var_name).sum(), np.ones(var_shape).sum())
            print('{} - get/set var: '.format(self.mod) + colored('PASS', 'green'))
            
            # check get_var / set_var at indices
            if self.bmi.grid._1d is not None:
                inds = self.bmi.grid._1d.inds
            else:
                inds = [0]
            self.bmi.set_value_at_indices(var_name, inds, np.ones_like(inds)*78) 
            self.assertEqual(self.bmi.get_value_at_indices(var_name, inds).sum(), np.ones_like(inds).sum()*78)
            print('{} - get/set var at indices: '.format(self.mod) + colored('PASS', 'green'))

            # finalize
            self.bmi.finalize()
            print('{} - finalize model: '.format(self.mod) + colored('PASS', 'green'))
        
        # cleanup
        finally:
            if self.outputDir: shutil.rmtree(self.outputDir)
            if self.mod != 'WFL':
                if self.bmi.initialized:
                    os.unlink(self.bmi._config_fn)
            if self.mod == 'CMF':
                if hasattr(self.bmi, '_mapdir'):
                    for fn in glob.glob(join(self.bmi._mapdir, '*tmp*')):
                        os.unlink(fn)
                    ext = ("*.txt","*.log","*.nc",)
                    for e in ext:
                        for fn in glob.glob(join(dirname(config_fn), e)):
                            os.unlink(fn)

        




if __name__ == "__main__":
    if len(sys.argv) > 1:
        MyTest.mod = sys.argv.pop()
    #     # MyTest.USERNAME = sys.argv.pop()
    unittest.main()


import unittest
from wflow.wflow_bmi import wflowbmi_csdms as wbmi
from pcrglobwb_bmi_v203.pcrglobwb_bmi import pcrglobwbBMI as pbmi
from glofrim import EBmi as abstact_bmi

class TestBmiObject(unittest.TestCase):

    def test_attrs(self):
        # bmi = wbmi()
        bmi = pbmi()
        for attr in abstact_bmi.__abstractmethods__:
            print(attr)
            self.assertTrue(hasattr(bmi, attr))

if __name__ == '__main__':
    unittest.main()
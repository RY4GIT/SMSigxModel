import unittest
import numpy as np
import pandas as pd
import sys

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/cfe/py_cfe")
import cfe
import bmi_cfe

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs")
from spotpy_cfe import spot_setup
from salib_cfe import SALib_CFE
from salib_cfe import salib_cfe_interface

class TestSALib(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('set up class')
        method_SALib = {
            'Morris': {'method': 'Morris', 'plot': 'EET', 'N': 1, 'n_levels': 2},
            # N=500, n_levels=4, total run = 8500 is ideal, N=3, n_levels=4 for a test run
            'Sobol': {'method': 'Sobol', 'plot': 'STS1', 'n': 2}
            # n=250 is ideal, n=2 for a test run

        }

        out_path= r'../4_out/sensitivity_analysis/Mahurangi/test'
        config_path_CFE= r'../2_data_input/unit_test/config_cfe.json'
        config_path_SALib= r'../2_data_input/unit_test/SALib_config.json'
        method_SALib=method_SALib['Morris']
        like_SALib='NashSutcliffe'
        var_measure_SALib='Flow'

        # Preperation
        cfe_instance = bmi_cfe.BMI_CFE(config_path_CFE)

        # Implementation
        cls._salib_experiment = SALib_CFE(
            cfe_instance=cfe_instance,
            config_path=config_path_SALib,
            method_SALib=method_SALib,
            like_measure=like_SALib,
            var_measure=var_measure_SALib,
        )

    @classmethod
    def tearDownClass(cls):
        print('tear down class')
        # cls.destroy()

    def setUp(self):
        print('setUp')
        self._salib_experiment.run()

    def tearDown(self):
        print('teadDown')

    def test_run_Sobol(self):
        print('test run(method=Sobol)')
        self.assertTrue(np.minimum(self._salib_experiment.Y) > 0.5)

    def test_run_Sobol(self):
        print('test run(method=Morris)')

if __name__ == '__main__':
    unittest.main()
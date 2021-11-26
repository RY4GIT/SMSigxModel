import os
import sys
sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs")
import numpy as np

import cfe
from spotpy_cfe import spot_setup
import spotpy
from salib_cfe import salib_cfe

from SALib.sample import morris as morris_s
from SALib.analyze import morris as morris_a
from SALib.test_functions import Ishigami

# specify current directory create output directory if it does not exist
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model")
os.getcwd()
out_file_path = '../4_out/test/'
if not os.path.exists(out_file_path):
    os.mkdir(out_file_path)
data_file_path = '../2_data_input/test'

def main(runtype):

    if runtype == "SALib":
        # test
        cfe1 = cfe.CFE(os.path.join(data_file_path, 'cat_58_config_cfe.json'))
        problem = {
            'num_vars':3,
            'names': ['maxsmc', 'wltsmc', 'init_soil_storage'],
            'bounds': [[0, 1.0],
                       [0, 1.0],
                       [0, 1.0]]
        }
        param_values = morris_s.sample(problem, 1)

        Y = np.zeros([param_values.shape[0]])
        for i, X in enumerate(param_values):
            Y[i] = salib_cfe(X, problem['names'], cfe1)

        Si = morris_a.analyze(problem, param_values, Y, print_to_console=True)
        print(Si['mu'])

    if runtype == "SPOTPy":
        # Initialize
        cfe1 = cfe.CFE(os.path.join(data_file_path, 'cat_58_config_cfe.json'))
        out_fn_sa = out_file_path + 'results'

        # Select number of maximum repetitions
        # Check out https://spotpy.readthedocs.io/en/latest/Sensitivity_analysis_with_FAST/ to determine an appropriate number of repetitions
        rep = 7

        # Start a sensitivity analysis
        sampler = spotpy.algorithms.fast(spot_setup(cfe_input=cfe1), dbname=out_fn_sa, dbformat='csv', save_sim = False)
        sampler.sample(rep)

        # Load the results gained with the fast sampler
        results = spotpy.analyser.load_csv_results(out_fn_sa)

        # Example plot to show the sensitivity index of each parameter
        spotpy.analyser.plot_fast_sensitivity(results, number_of_sensitiv_pars=2)

        # Example to get the sensitivity index of each parameter
        SI = spotpy.analyser.get_sensitivity_of_fast(results)

if __name__ == '__main__':
    main(runtype = "SALib")



"""
cfe1.initialize()
cfe1.update()
cfe1.update_until(4)
cfe1.finalize()
"""
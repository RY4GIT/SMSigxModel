# A main module to run various analysis with CFE model

import os
import sys
import numpy as np
import pandas as pd
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import spotpy
import cProfile, pstats, io
from pstats import SortKey
# https://docs.python.org/3/library/profile.html#module-cProfile

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/cfe/py_cfe")
import cfe
import bmi_cfe

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/SMSig")
from sig_seasontrans import SMSig

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs")
from spotpy_cfe import spot_setup
from salib_cfe import SALib_CFE

# sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/glue")
from glue_cfe import MyGLUE

# specify current directory create output directory if it does not exist
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model")
os.getcwd()
out_path = '../4_out/Mahurangi/'
if not os.path.exists(out_path):
    os.mkdir(out_path)
data_file_path = '../2_data_input/Mahurangi'

# from numba import jit
# @jit

def main(runtype):

    if runtype == "NOAA_CFE":
        cfe_instance = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'config_cfe.json'))
        cfe_instance.initialize()
        cfe_instance.run_unit_test(plot=True, print_fluxes=False)
        cfe_instance.finalize()

    if runtype == "SALib":
        print('Start sensitivity analysis')

        # preparation & sampling parameters
        cfe_instance = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'full', 'config_cfe.json'))
        problem = {
            'num_vars':6,
            'names': ['smcmax', 'wltsmc', 'K_lf', 'K_nash', 'D', 'bb'],
            'bounds': [[0, 1.0],
                       [0, 1.0],
                       [0.000001, 100],
                       [0.00001,30],
                       [1,10],
                       [1,6]]
        }

        salib_experiment = SALib_CFE(
            cfe_instance=cfe_instance, problem=problem, SAmethod='Sobol', out_path=out_path
        )
        salib_experiment.run()
        salib_experiment.plot(plot_type='STS1')


    if runtype == "SPOTPy":
        # Initialize
        cfe1 = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'config_cfe.json'))
        cfe1.initialize()
        out_fn_sa = out_path + 'results'

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

    if runtype == "GLUE":
        # Initialize
        cfe1 = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'full', 'config_cfe.json'))
        cfe1.initialize()
        out_fn_sa = out_path + 'results'

        # Start GLUE
        nrun = 3
        glue1 = MyGLUE(cfe_input = cfe1, out_path=out_path, nrun=nrun)
        glue1.simulation()
        glue1.post_process()
        glue1.to_csv()
        glue1.plot()

    if runtype == "Seasonsig":

        # Get the comparison data
        cfe1 = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'full', 'config_cfe.json'))
        cfe1.initialize()
        obs0 = cfe1.load_unit_test_data()
        obs = obs0[["Time", "Soil Moisture Content"]]
        obs["Time"] = pd.to_datetime(obs["Time"], format="%d-%b-%Y %H:%M:%S")  # Works specifically for Mahurangi data

        sig = SMSig(ts_time=obs["Time"].to_numpy(), ts_value=obs["Soil Moisture Content"].to_numpy())
        sig.movmean()
        sig.detrend()
        t_valley = sig.calc_sinecurve()
        season_trans1 = sig.calc_seasontrans(t_valley=t_valley)

    """
    if runtype == "cfe_CUAHSI":
        cfe1 = cfe.CFE(os.path.join(data_file_path, 'cat_58_config_cfe.json'))
        cfe1.initialize()
        cfe1.update()
        cfe1.update_until(4)
        cfe1.finalize()
    """

if __name__ == '__main__':

    measuretime = True
    # measure the time
    if measuretime:
        pr = cProfile.Profile()
        pr.enable()

    main(runtype = "GLUE")

    # measure the time
    if measuretime:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        ps.dump_stats('runtime.txt')

# snakeviz "G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\3_model\runtime.txt"
# visualize to type the above in terminal
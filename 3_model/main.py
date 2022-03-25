# A main module to run various analysis with CFE model

import os
import sys
import numpy as np
import pandas as pd
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
import json

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
out_path = '..\\4_out\\Mahurangi'
if not os.path.exists(out_path):
    os.mkdir(out_path)
data_file_path = '..\\2_data_input\\Mahurangi\\full'

# from numba import jit
# @jit

def main(runtype, nrun):

    if runtype == "NOAA_CFE":

        input_config = {
            "forcing_file": "..\\2_data_input\\Mahurangi\\full\\mahurangi_1998_2001.csv",
            "catchment_area_km2": 46.65,
            "soil_params": {
                "bb": 14.658880935233976,
                "mult": 1000.0,
                "satdk": 0.013447137666802389,
                "satpsi": 0.09187291127012216,
                "slop": 0.0002740918314719565,
                "smcmax": 0.6757270053046729,
                "wltsmc": 0.3689068227673663, "D": 0.87,
                "exponent_primary": 1.0,
                "coeff_secondary": 1.0174974070860094,
                "exponent_secondary": 1.1833082840205675
            },
            "max_gw_storage": 221.96089432793116,
            "Cgw": 0.9556998526934358,
            "expon": 1.073458751801099,
            "K_lf": 46.875000531249995,
            "K_nash": 0.22440411667008253,
            "nash_storage": [0.0, 0.0],
            "giuh_ordinates": [0.1, 0.35, 0.2, 0.14, 0.1, 0.06, 0.05],
            "refkdt": 2.325787803305572,
            "trigger_z_m": 0.8340054706423846,
            "fc_atm_press_fraction": 0.33,
            "stand_alone": 1,
            "unit_test": 1,
            "compare_results_file": "G:\\Shared drives\\Ryoko and Hilary\\SMSigxModel\\analysis\\2_data_input\\Mahurangi\\full\\test_sm_basinavg.csv",
        }

        input_path = os.path.join(data_file_path, 'full', 'config_cfe.json')
        with open(input_path, 'w') as outfile:
            json.dump(input_config, outfile)

        with open(input_path) as outfile:
            loaded_data = json.load(outfile)

        cfe_instance = bmi_cfe.BMI_CFE(input_path)
        cfe_instance.initialize()
        cfe_instance.run_unit_test(plot=True, print_fluxes=False)
        cfe_instance.finalize()

    if runtype == "SALib":
        print('Start sensitivity analysis')

        # preparation & sampling parameters
        cfe_instance = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'full', 'config_cfe.json'))
        problem = {
            'num_vars':16,
            'names': ['bb',
                      'satdk',
                      'satpsi',
                      'slop',
                      'smcmax',
                      'wltsmc',
                      'D',
                      'coeff_secondary',
                      'exponent_secondary',
                      'max_gw_storage',
                      'Cgw',
                      'expon',
                      'K_nash',
                      'refkdt',
                      'trigger_z_m',
                      'fc_atm_press_fraction'
                      ],
            'bounds': [[2, 15], # bb
                       [0, 1],  # satdk
                       [0.02, 0.78], # satpsi
                       [0,1], # slop
                       [0.33, 0.7], # smcmax
                       [0,0.5], #wltsmc
                       [0.01, 2], # D
                       [0.01, 3], #coeff_secondary
                       [1,8], #exponent_secondary
                       [10,250], # max_gw_storage
                       [0.01,3], # Cgw
                       [1,8], # expon
                       [0,1], #K_nash
                       [0.1,4], #refkdt
                       [0.01,0.87], # trigger_z_m
                       [0.01,0.33] #fc_atm_press_fraction
                       ]
        }

        """
        salib_experiment = SALib_CFE(
            cfe_instance=cfe_instance, problem=problem, SAmethod='Sobol', out_path=out_path
        )
        salib_experiment.run()
        salib_experiment.plot(plot_type='STS1')
        """

        out_path_salibexp = '..\\4_out\\Mahurangi\\salidexp_id1'
        if not os.path.exists(out_path_salibexp):
            os.mkdir(out_path_salibexp)

        salib_experiment = SALib_CFE(
            cfe_instance=cfe_instance, problem=problem, SAmethod='Morris', out_path=out_path_salibexp
        )
        salib_experiment.run()
        salib_experiment.plot(plot_type='EET')

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
        cfe1 = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'config_cfe.json'))
        cfe1.initialize()

        out_path_glueexp = '..\\4_out\\Mahurangi\\'
        if not os.path.exists(out_path_glueexp):
            os.mkdir(out_path_glueexp)

        # Start GLUE
        glue1 = MyGLUE(cfe_input = cfe1, out_path=out_path_glueexp, nrun=nrun, calib_case=3)
        glue1.simulation()
        glue1.post_process()
        glue1.to_csv()
        # glue1.plot(plot_type="dotty")
        glue1.plot(plot_type="dotty_interaction")
        # glue1.plot(plot_type="param_hist")
        # glue1.plot(plot_type="timeseries")

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

    # o = open(os.path.join(out_path, 'log.txt'), 'w')
    main(runtype = "GLUE")
    #o.close()

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
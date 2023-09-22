# A main module to run various analysis with CFE model

# Import libraries
import os
import sys
import numpy as np
import pandas as pd
import sys

if not sys.warnoptions:
    import warnings

    warnings.simplefilter("ignore")
import json
import multiprocessing as mp
import snakeviz
import spotpy
import cProfile, pstats, io
from pstats import SortKey

# https://docs.python.org/3/library/profile.html#module-cProfile

sys.path.append(
    "G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/cfe/cfe_py"
)
import cfe
import bmi_cfe

sys.path.append(
    "G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/SMSig"
)
from sig_seasontrans import SMSig

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs")
from spotcfe_py import spot_setup
from salib_cfe import SALib_CFE

# sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/glue")
from glue_cfe import MyGLUE

# Specify current directory create output directory if it does not exist
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model")
os.getcwd()


def SALib_excel_to_config(config_excel_path_SALib, config_json_path_SALib):
    # This function reads configuration setting for SALib sensitivity analysis
    # from Excel file to json file format

    # Read the data
    df = pd.read_excel(config_excel_path_SALib)
    df_param_to_calibrate = df[df["calibrate"] == 1]

    # Convert it to dictionary
    config_SALib = dict()
    config_SALib["num_vars"] = df_param_to_calibrate.shape[0]
    config_SALib["names"] = [None] * df_param_to_calibrate.shape[0]
    config_SALib["bounds"] = [[np.nan, np.nan]] * df_param_to_calibrate.shape[0]
    for i in df_param_to_calibrate.index:
        config_SALib["names"][i] = df_param_to_calibrate["name"][i]
        config_SALib["bounds"][i] = [
            df_param_to_calibrate["lower_bound"][i],
            df_param_to_calibrate["upper_bound"][i],
        ]

        if (
            df_param_to_calibrate["upper_bound"][i]
            <= df_param_to_calibrate["lower_bound"][i]
        ):
            warnings.warn(
                "The upper bound is smaller than lower bound for the parameter %s"
                % {config_SALib["names"][i]}
            )

    # Output the data
    with open(config_json_path_SALib, "w") as outfile:
        json.dump(config_SALib, outfile, indent=4)


# Create an instance
def main(
    runtype,
    out_path="../4_out/",
    config_path_CFE="",
    config_path_SALib="",
    method_SALib=None,
    like_SALib="",
    var_measure_SALib="",
    config_path_GLUE="",
    nrun_GLUE=1,
    glue_calib_case=1,
):
    if runtype == "run_CFE":
        # To simply run the CFE model
        print("Run the CFE model")

        cfe_instance = bmi_cfe.BMI_CFE(config_path_CFE)
        cfe_instance.initialize()
        cfe_instance.run_unit_test(plot=True, print_fluxes=False)
        cfe_instance.finalize(print_mass_balance=True)

    if runtype == "SALib":
        # To implement sensitivity analysis with SALib. Currently this module supports Morris and Sobol analysis
        print("Start sensitivity analysis with ** %s **" % method_SALib["method"])

        # Preperation
        cfe_instance = bmi_cfe.BMI_CFE(config_path_CFE)
        if not os.path.exists(out_path):
            os.mkdir(out_path)

        # Implementation
        salib_experiment = SALib_CFE(
            cfe_instance=cfe_instance,
            config_path=config_path_SALib,
            method_SALib=method_SALib,
            like_measure=like_SALib,
            var_measure=var_measure_SALib,
            out_path=out_path,
        )
        salib_experiment.run()
        salib_experiment.plot()

    if runtype == "SPOTPy":
        # To implement sensitivity analysis with SPOTPy.
        print("Start sensitivity analysis")

        # Initialize
        cfe1 = bmi_cfe.BMI_CFE(config_path_CFE)
        cfe1.initialize()
        out_fn_sa = out_path + "results"

        # Select number of maximum repetitions
        # Check out https://spotpy.readthedocs.io/en/latest/Sensitivity_analysis_with_FAST/ to determine an appropriate number of repetitions
        rep = 7

        # Start a sensitivity analysis
        sampler = spotpy.algorithms.fast(
            spot_setup(cfe_input=cfe1), dbname=out_fn_sa, dbformat="csv", save_sim=False
        )
        sampler.sample(rep)

        # Load the results gained with the fast sampler
        results = spotpy.analyser.load_csv_results(out_fn_sa)

        # Example plot to show the sensitivity index of each parameter
        spotpy.analyser.plot_fast_sensitivity(results, number_of_sensitiv_pars=2)

        # Example to get the sensitivity index of each parameter
        SI = spotpy.analyser.get_sensitivity_of_fast(results)

    if runtype == "GLUE":
        # To implement sensitivity analysis with GLUE.
        print("Start GLUE analysis")

        # Initialize
        cfe_instance = bmi_cfe.BMI_CFE(config_path_CFE)
        cfe_instance.initialize()

        if not os.path.exists(out_path):
            os.mkdir(out_path)

        # Start GLUE
        glue_instance = MyGLUE(
            cfe_input=cfe_instance,
            out_path=out_path,
            config_path=config_path_GLUE,
            nrun=nrun_GLUE,
            calib_case=glue_calib_case,
        )
        glue_instance.simulation()
        glue_instance.post_process()

        # Output the results
        glue_instance.to_csv()
        glue_instance.plot(plot_type="dotty")
        glue_instance.plot(plot_type="dotty_interaction")
        glue_instance.plot(plot_type="param_hist")
        glue_instance.plot(plot_type="timeseries")

    if runtype == "Seasonsig":
        # To test soil moisture signature
        print("Start seasonal signature test")

        # Get the comparison data
        cfe1 = bmi_cfe.BMI_CFE(os.path.join(data_file_path, "full", "config_cfe.json"))
        cfe1.initialize()
        obs0 = cfe1.load_unit_test_data()
        obs = obs0[["Time", "Soil Moisture Content"]]
        obs["Time"] = pd.to_datetime(
            obs["Time"], format="%d-%b-%Y %H:%M:%S"
        )  # Works specifically for Mahurangi data

        sig = SMSig(
            ts_time=obs["Time"].to_numpy(),
            ts_value=obs["Soil Moisture Content"].to_numpy(),
        )
        sig.movmean()
        sig.detrend()
        t_valley = sig.calc_sinecurve()
        season_trans1 = sig.calc_seasontrans(t_valley=t_valley)


if __name__ == "__main__":
    # measure the running time
    measuretime = True
    if measuretime:
        pr = cProfile.Profile()
        pr.enable()

    # ===============================================
    # =========== RUN CFE ==============
    # ===============================================

    # main(
    #     runtype="run_CFE",
    #     out_path='../4_out/unit_test/',
    #     config_path_CFE='../2_data_input/unit_test/config_cfe.json',
    #     )

    # ===============================================
    # =========== SENSITIVITY ANALYSIS ==============
    # ===============================================

    # To implement sensitivity analysis with SALib
    # Specify the sensitivity analysis method, plotting method, and number of runs here
    method_SALib = {
        "Morris": {
            "method": "Morris",
            "plot": "EET",
            "N": 500,
            "n_levels": 4,
        },  # n=250 is ideal, n=2 for a test run
        "Sobol": {
            "method": "Sobol",
            "plot": "STS1",
            "n": 250,
        },  # N=500, n_levels=4, total run = 8500 is ideal, N=3, n_levels=4 for a test run
        "stability_test": {"method": "stability_test", "plot": "STS1", "n": 250},
    }

    SALib_excel_to_config(
        config_excel_path_SALib="../2_data_input/unit_test/SALib_config_wide.xlsx",
        config_json_path_SALib="../2_data_input/unit_test/SALib_config_wide.json",
    )

    # Run the analysis
    main(
        runtype="SALib",
        out_path="../4_out/sensitivity_analysis/Mahurangi/test",
        config_path_CFE="../2_data_input/unit_test/short_config_cfe.json",
        config_path_SALib="../2_data_input/unit_test/SALib_config_wide.json",
        method_SALib=method_SALib["stability_test"],
        like_SALib="NashSutcliffe",
        var_measure_SALib="Flow",  # 'Flow' (discharge in meter), 'Soil Moisture Content' soil moisture content in fraction
    )

    # ===============================================
    # =========== SENSITIVITY ANALYSIS ==============
    # ===============================================
    main(
        runtype="GLUE",
        out_path="../4_out/sensitivity_analysis/Mahurangi/test/GLUE",
        config_path_CFE="../2_data_input/unit_test/short_config_cfe.json",
        config_path_GLUE="../2_data_input/unit_test/SALib_config_wide.json",
        nrun_GLUE=1,
    )

    # main(
    #     runtype="SALib",
    #     out_path='../4_out/sensitivity_analysis/Mahurangi/test',
    #     config_path_CFE='../2_data_input/unit_test/config_cfe.json',
    #     config_path_SALib='../2_data_input/unit_test/SALib_config.json',
    #     method_SALib=method_SALib['Sobol'],
    #     like_SALib = 'NashSutcliffe',
    #     var_measure_SALib='Flow'
    # )

    # ===============================================
    # =========== GLUE ANALYSIS ==============
    # ===============================================

    # main(runtype="GLUE", nrun=5000, glue_calib_case=1, out_path='..\\4_out\\Mahurangi\\',cfe_json_fn='config_cfe.json')
    # data_file_path = '..\\2_data_input\\Mahurangi\\full'
    # main(runtype="SALib", out_path='..\\4_out\\Mahurangi\\SALib_id1', data_file_path=data_file_path)
    # nrun = 10000
    # main(runtype = "GLUE", nrun=nrun, glue_calib_case=1, out_path= '..\\4_out\\Mahurangi\\exp_id9') #NSE_Q
    # main(runtype="GLUE", nrun=nrun, glue_calib_case=2, out_path= '..\\4_out\\Mahurangi\\exp_id10') #KGE_Q
    # main(runtype="GLUE", nrun=nrun, glue_calib_case=5, out_path= '..\\4_out\\Mahurangi\\exp_id11') #sesasonSM
    # main(runtype="GLUE", nrun=nrun, glue_calib_case=3, out_path= '..\\4_out\\Mahurangi\\exp_id12') #KGE_SM
    # main(runtype="GLUE", nrun=nrun, glue_calib_case=6, out_path= '..\\4_out\\Mahurangi\\exp_id6') #multi

    """
    process1 = mp.Process(target=main, kwargs={'runtype': "GLUE", 'nrun': nrun, 'glue_calib_case': 1,
                                               'out_path': '..\\4_out\\Mahurangi\\exp_id9', 'cfe_json_fn': 'config_cfe_core1.json'})
    process2 = mp.Process(target=main, kwargs={'runtype': "GLUE", 'nrun': nrun, 'glue_calib_case': 2,
                                               'out_path': '..\\4_out\\Mahurangi\\exp_id10', 'cfe_json_fn': 'config_cfe_core2.json'})
    process3 = mp.Process(target=main, kwargs={'runtype': "GLUE", 'nrun': nrun, 'glue_calib_case': 5,
                                               'out_path': '..\\4_out\\Mahurangi\\exp_id11', 'cfe_json_fn': 'config_cfe_core3.json'})
    process4 = mp.Process(target=main, kwargs={'runtype': "GLUE", 'nrun': nrun, 'glue_calib_case': 3,
                                               'out_path': '..\\4_out\\Mahurangi\\exp_id12', 'cfe_json_fn': 'config_cfe_core4.json'})
    # process2 = mp.Process(target=main)

    process1.start()
    process2.start()
    process3.start()
    process4.start()
    """

    # measure the time
    if measuretime:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        # print(s.getvalue())
        ps.dump_stats("runtime.txt")

    # "C:\Users\flipl\anaconda3\envs\CFE\Scripts\snakeviz.exe"
    # snakeviz "G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\3_model\runtime.txt"
    # visualize to type the above in terminal

# A main module to run various analysis with CFE model

# Import libraries
import os
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import spotpy

import sys
sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/libs/cfe/py_cfe")
import cfe
import bmi_cfe

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/libs/SMSig")
from sig_seasontrans import SMSig

# Specify current directory create output directory if it does not exist
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/5_GLUE_model")
os.getcwd()

def main(out_path='', config_path_CFE='', config_path_GLUE='', eval_criteria=dict()):

    in_path = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\6_out\Mahurangi\ex1\paramter_priori.csv"
    config_temp = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\2_data_input\Mahurangi\parameters\config_cfe_0.json"
    out_path = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\6_out\seasonsig_test"
    i_run = 9

    # Read parameters from excel sheet and create a config file
    config_all_runs = pd.read_csv(in_path)
    config_target_runs = config_all_runs.iloc[i_run]

    with open(config_temp, 'r') as outfile:
        config_temp = json.load(outfile)

    for i in range(len(config_target_runs)):
        if config_target_runs.index[i] in ['bb', 'satdk', 'slop', 'satpsi', 'smcmax', 'wltsmc', 'D']:
            config_temp["soil_params"][config_target_runs.index[i]] = config_target_runs[i]
        else:
            config_temp[config_target_runs.index[i]] = config_target_runs[i]

    with open(os.path.join(out_path, 'config_CFE.json'), 'w') as out_file:
        json.dump(config_temp, out_file)

    # Run the CFE based on the config file
    cfe_instance = bmi_cfe.BMI_CFE(os.path.join(out_path, 'config_CFE.json'))
    cfe_instance.initialize()
    sim0 = cfe_instance.run_unit_test(plot=False, warm_up=True)
    obs0 = cfe_instance.load_unit_test_data()

    # Get the results
    var_name = 'Soil Moisture Content'
    sim = sim0[["Time", var_name]].copy()
    sim["Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S")  # Works specifically for CFE

    # Get the comparison data
    obs = obs0[["Time", var_name]].copy()
    obs["Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M")  # Works specifically for Mahurangi data
    # obs["Time"] = pd.to_datetime(obs["Time"], format="%d-%m-%Y %H:%M:%S")

    # Merge observed and simulated timeseries
    df = pd.merge_asof(sim, obs, on="Time")

    sim_synced = pd.DataFrame()
    obs_synced = pd.DataFrame()
    sim_synced[var_name] = df[var_name + "_x"].copy()
    obs_synced[var_name] = df[var_name + "_y"].copy()

    # Evaluate using seasonal soil moisture signature
    sig_obs = SMSig(
        ts_time=df["Time"].to_numpy(),
        ts_value=obs_synced[var_name].to_numpy(),
        plot_results=True,
        plot_label="obs"
    )
    # sig_obs.detrend() # TODO:debug
    sig_obs.movmean()
    t_valley = sig_obs.calc_sinecurve()
    season_trans_obs, start_dates_obs, end_dates_obs = sig_obs.calc_seasontrans(t_valley=t_valley)

    sig_sim = SMSig(
        ts_time=df["Time"].to_numpy(),
        ts_value=sim_synced[var_name].to_numpy(),
        plot_results=True,
        plot_label="sim"
    )
    sig_sim.movmean()
    season_trans_sim, start_dates_sim, end_dates_sim = sig_sim.calc_seasontrans(t_valley=t_valley)

    # Get the deviations in seasonal transition dates between simulated and observed timeseries
    diff = season_trans_sim - season_trans_obs
    metric_value = abs(np.nanmean(diff, axis=0))
    print(diff)
    print(metric_value)

    # Plot out the results
    df_obs = obs_synced
    df_sim = sim_synced
    obs_label = 'Observed'
    sim_label = 'Simulated'
    obs_color = '#1f77b4'
    sim_color = '#ff7f0e'
    y_label = 'Volumetric Soil Moisture Content [m^3/m^3]'
    title = 'Soil moisture and seasonal transition signatures'
    fn = 'timeseries.png'

    f2 = plt.figure(figsize=(15, 5))
    ax2 = f2.add_subplot()
    ax2.plot(sig_obs.tt.index, sig_obs.tt.values, alpha=1, label=obs_label, color=obs_color)
    ax2.plot(sig_sim.tt.index, sig_sim.tt.values, alpha=1, label=sim_label, color=sim_color)
    # ax2.plot(df["Time"].values, df_obs[var_name].values, alpha=1, label=obs_label, color=obs_color)
    # ax2.plot(df["Time"].values, df_sim[var_name].values, alpha=1, label=sim_label, color=sim_color)
    for i in range(len(start_dates_obs)):
        ax2.axvline(x=start_dates_obs[i], color=obs_color, label=None, alpha=0.5, linestyle='-')
    for i in range(len(end_dates_obs)):
        ax2.axvline(x=end_dates_obs[i], color=obs_color, label=None, alpha=0.5, linestyle='--')
    for i in range(len(start_dates_sim)):
        ax2.axvline(x=start_dates_sim[i], color=sim_color, label=None, alpha=0.5, linestyle='-')
    for i in range(len(end_dates_sim)):
        ax2.axvline(x=end_dates_sim[i], color=sim_color, label=None, alpha=0.5, linestyle='--')
    ax2.set_xlabel('Time')
    ax2.set_ylabel(y_label)
    ax2.set_title(title)
    ax2.legend()
    f2.savefig(os.path.join(out_path, fn), dpi=600)


    # Save the results


if __name__ == '__main__':
    main()
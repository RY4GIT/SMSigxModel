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
import seaborn as sns
from math import log10

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
    config_temp = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\2_data_input\Mahurangi\parameters\ex1_config_cfe.json"
    out_path = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\6_out\post_eval_test"
    i_run = 3

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
    var_name = 'Flow'
    sim = sim0[["Time", var_name]].copy()
    sim["Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S")  # Works specifically for CFE

    # Get the comparison data
    obs = obs0[["Time", var_name]].copy()
    obs["Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M")  # Works specifically for Mahurangi data
    # obs["Time"] = pd.to_datetime(obs["Time"], format="%d-%m-%Y %H:%M:%S")

    obs_P = obs0[["Time", "Rainfall"]].copy()
    obs_P["Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M")

    # Merge observed and simulated timeseries
    df = pd.merge_asof(sim, obs, on="Time")
    df2 = pd.merge_asof(df, obs_P, on="Time")

    sim_synced = pd.DataFrame()
    obs_synced = pd.DataFrame()

    sim_synced[var_name] = df[var_name + "_x"].copy()
    obs_synced[var_name] = df[var_name + "_y"].copy()
    obs_synced["Rainfall"] = df2["Rainfall"].copy()

    print('stop')
    df_obs = obs_synced
    df_sim = sim_synced

    df_obs.set_axis(df["Time"], axis=0, inplace=True)
    df_sim.set_axis(df["Time"], axis=0, inplace=True)
    # Test some new functions

    # Run off ratio
    print(f'Runoff ratio (over the simulated period): {df_obs["Flow"].sum()/df_obs["Rainfall"].sum()}')
    print(f'Sum of the rainfall (m): {df_obs["Rainfall"].sum()}')
    print(f'Sum of the flow (m): {df_obs["Flow"].sum()}')

    df_obs_monthly = df_obs.resample('M').sum().copy()
    df_sim_monthly = df_sim.resample('M').sum().copy()
    df_obs_monthly.drop(df_obs_monthly.tail(1).index, inplace=True)
    df_obs_monthly.drop(df_obs_monthly.head(1).index, inplace=True)
    df_sim_monthly.drop(df_sim_monthly.tail(1).index, inplace=True)
    df_sim_monthly.drop(df_sim_monthly.head(1).index, inplace=True)

    df_obs_runoff_ratio = df_obs_monthly["Flow"]/ df_obs_monthly["Rainfall"]
    df_sim_runoff_ratio = df_sim_monthly["Flow"] / df_obs_monthly["Rainfall"]
    bias_runoff_ratio = df_sim_runoff_ratio-df_obs_runoff_ratio

    ## Variabiltiy index
    #  https://sebastiangnann.github.io/TOSSH_development/matlab/TOSSH_code/TOSSH_development/TOSSH_code/signature_functions/sig_VariabilityIndex.html
    df_obs_sorted = df_obs.sort_values(by=['Flow'], ascending=False).copy()
    df_sim_sorted = df_sim.sort_values(by=['Flow'], ascending=False).copy()
    percs = np.arange(10, 100, 10)
    indices_percs = np.round(len(df_obs_sorted)*percs/100);
    flow_values = df_obs_sorted['Flow'].values
    flow_percs = flow_values[indices_percs.astype(int)]
    recs = flow_percs > 0
    VariabilityIndex = np.std(np.log10(flow_percs[recs]))


    # Plot runoff ratio results
    df_results = pd.concat([df_obs_runoff_ratio, df_sim_runoff_ratio, bias_runoff_ratio], axis=1)
    new_name = ['observed_runoff_ratio', 'simulated_runoff_ratio', 'bias_runoff_ratio']
    for i in range(len(list(df_results))):
        df_results.rename(columns={list(df_results)[i]: new_name[i]}, inplace=True)

    df_results.to_csv(os.path.join(out_path, 'runoff_ratio.csv'), sep=',', index=True,
                                         encoding='utf-8', na_rep='nan')
    # Plot out the results
    obs_label = 'Observed'
    sim_label = 'Simulated'
    obs_color = '#1f77b4'
    sim_color = '#ff7f0e'

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    ax1 = fig.add_subplot(2, 1, 1)
    line1, = ax1.plot(df_obs_runoff_ratio, label=obs_label, color=obs_color)
    line2, = ax1.plot(df_sim_runoff_ratio, label=sim_label, color=sim_color)
    ax2 = fig.add_subplot(2, 1, 2)
    bp = sns.boxplot(data=df_results.transpose(), ax=ax2)
    fig.legend()
    xax = ax1.xaxis
    ax1.set_title("Results from CFE")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Monthly runoff ratio [Q/P]")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Bias in runoff ratio\n[Simulated - observed]")
    fig.autofmt_xdate()
    fig.savefig(os.path.join(out_path, 'runoff_ratio.png'))
    del fig, ax1, ax2


    # Save the results


if __name__ == '__main__':
    main()
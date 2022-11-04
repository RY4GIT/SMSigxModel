import sys
sys.path.append("../libs/cfe/py_cfe")
import cfe
import bmi_cfe
sys.path.append("../libs/SMSig")
from sig_seasontrans import SMSig

import spotpy
import os
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import shutil
from math import log10
from statistics import median

# Global function


# GLUE object
class MyGLUEPost(object):
    def __init__(self, out_path='./', config_path_CFE=''):
        print("Post evaluation of GLUE analysis")

        eval_criteria = {
        0: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'KGE', 'threshold': -1000},
        1: {'variable_to_analyze': 'Flow', 'metric': 'KGE', 'threshold':  -1000},
        2: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'season_transition', 'threshold': -1000}
        }

        ## Variabiltiy index
        #  https://sebastiangnann.github.io/TOSSH_development/matlab/TOSSH_code/TOSSH_development/TOSSH_code/signature_functions/sig_VariabilityIndex.html

        eval_criteria_monthly = {
            0: {'variable_to_analyze': 'Flow', 'metric': 'Q_mean'},
            1: {'variable_to_analyze': 'Flow', 'metric': 'high_flow_freq'},
            2: {'variable_to_analyze': 'Flow', 'metric': 'RR'}
        }

        self.out_path = out_path  # Output folder path
        self.var_names = ["Flow", "Soil Moisture Content"]  # Variables to be analyzed
        self.config_path_CFE = config_path_CFE

        # Evaluation criteria (multi-criteria allowed)
        self.eval_criteria = eval_criteria
        # print(f"A number of criterion: {len(eval_criteria)}")
        self.eval_names = []
        for i in range(len(eval_criteria)):
            if eval_criteria[i]["metric"] == 'season_transition':
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (d2w_start)')
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (d2w_end)')
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (w2d_start)')
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (w2d_end)')
            else:
                self.eval_names.append(f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]}')
            # print(f'[{i + 1}] {eval_criteria[i]["metric"]}-based analysis on {eval_criteria[i]["variable_to_analyze"]}')

        self.eval_criteria_monthly = eval_criteria_monthly

    def simulation(self, behavioral_param):

        # ===============================================================
        # The main model run
        # ===============================================================

        # Preparations
        nth_run = behavioral_param[0]
        sampled_params_set = behavioral_param[1]
        print(f'Processing {nth_run}')

        # ===============================================================
        # Write the randomly-generated parameters to the config json file
        # ===============================================================

        # CFE model instance
        template_config_CFE = self.config_path_CFE
        target_config_CFE = os.path.join(r"..\2_data_input\Mahurangi\parameters",
                                         f"config_cfe_{nth_run}.json")
        shutil.copyfile(template_config_CFE, target_config_CFE)

        # Get the model config file
        with open(target_config_CFE) as data_file:
            self.cfe_cfg = json.load(data_file)

        # Overwrite the model config file
        for i, target_para in enumerate(sampled_params_set.index):
            if target_para in ['bb', 'satdk', 'slop', 'satpsi', 'smcmax', 'wltsmc', 'D']:
                self.cfe_cfg["soil_params"][target_para] = sampled_params_set[target_para]
            else:
                self.cfe_cfg[target_para] = sampled_params_set[target_para]

        # Save the config file with new parameters
        with open(target_config_CFE, 'w') as out_file:
            json.dump(self.cfe_cfg, out_file)

        # ===============================================================
        # Actual model run
        # ===============================================================
        self.myCFE = bmi_cfe.BMI_CFE(target_config_CFE)
        self.myCFE.initialize()
        sim0 = self.myCFE.run_unit_test(plot=False, warm_up=True)
        obs0 = self.myCFE.load_unit_test_data()

        # ===============================================================
        # Retrieve the simulated and observed data
        # ===============================================================
        sim_synced = pd.DataFrame()
        obs_synced = pd.DataFrame()

        # Get the results
        for var_name in self.var_names:
            # Get the simulated data
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

            sim_synced[var_name] = df[var_name + "_x"].copy()
            obs_synced[var_name] = df[var_name + "_y"].copy()
            obs_synced["Rainfall"] = df2["Rainfall"].copy()


        # ===============================================================
        # Evalute the outputs (whole-period metrics)
        # ===============================================================

        # Preparation
        eval_result_for_a_run = []
        behavioral_flag = [False] * len(self.eval_criteria)

        # Loop for all evaluation metrics (multi-criteria).
        # Calculate the metrics and judge behavioral vs. non-behavioral
        for i in range(len(self.eval_criteria)):

            # Nash-Sutcliffe scores
            if self.eval_criteria[i]['metric'] == "NSE":
                metric_value = spotpy.objectivefunctions.nashsutcliffe(
                    obs_synced[self.eval_criteria[i]['variable_to_analyze']],
                    sim_synced[self.eval_criteria[i]['variable_to_analyze']]
                )
                if metric_value > self.eval_criteria[i]['threshold']:
                    behavioral_flag[i] = True

            # Kling-Gupta Efficiency scores
            # 　Kling-Gupta efficiencies range from -Inf to 1. Essentially, the closer to 1, the more accurate the model is
            elif self.eval_criteria[i]['metric'] == "KGE":
                metric_value = spotpy.objectivefunctions.kge(
                    obs_synced[self.eval_criteria[i]['variable_to_analyze']],
                    sim_synced[self.eval_criteria[i]['variable_to_analyze']]
                )
                if metric_value > self.eval_criteria[i]['threshold']:
                    behavioral_flag[i] = True

            # Seasonal transition dates
            elif self.eval_criteria[i]['metric'] == "season_transition":
                # Calculate metrics for OBSERVED timeseries as a baseline performance. Run only once
                # if n == 0:
                sig_obs = SMSig(
                    ts_time=df["Time"].to_numpy(),
                    ts_value=obs_synced[self.eval_criteria[i]['variable_to_analyze']].to_numpy(),
                    plot_results=False,
                    plot_label="obs"
                )
                # sig_obs.detrend() # TODO:debug
                sig_obs.movmean()
                t_valley = sig_obs.calc_sinecurve()
                season_trans_obs, _, _ = sig_obs.calc_seasontrans(t_valley=t_valley)

                # Calculate metrics for SIMULATED timeseries
                sig_sim = SMSig(
                    ts_time=df["Time"].to_numpy(),
                    ts_value=sim_synced[self.eval_criteria[i]['variable_to_analyze']].to_numpy(),
                    plot_results=False,
                    plot_label="sim"
                )
                sig_sim.movmean()
                season_trans_sim, _, _ = sig_sim.calc_seasontrans(t_valley=t_valley)

                # Get the deviations in seasonal transition dates between simulated and observed timeseries
                diff = season_trans_sim - season_trans_obs
                metric_value = abs(np.nanmean(diff, axis=0))
                if all(metric_value < self.eval_criteria[i]['threshold']):
                    behavioral_flag[i] = True

            # Store evaluation metrics for all criteria for one run
            eval_result_for_a_run.append(metric_value)

        # ===============================================================
        # Evalute the outputs (monthly metrics)
        # ===============================================================

        df_obs = obs_synced
        df_sim = sim_synced

        df_obs.set_axis(df["Time"], axis=0, inplace=True)
        df_sim.set_axis(df["Time"], axis=0, inplace=True)

        df_obs_monthly = df_obs.resample('M').sum().copy()
        df_sim_monthly = df_sim.resample('M').sum().copy()
        df_obs_monthly.drop(df_obs_monthly.tail(1).index, inplace=True)
        df_obs_monthly.drop(df_obs_monthly.head(1).index, inplace=True)
        df_sim_monthly.drop(df_sim_monthly.tail(1).index, inplace=True)
        df_sim_monthly.drop(df_sim_monthly.head(1).index, inplace=True)

        median_flow_obs = median(df_obs['Flow'])
        high_flow_obs = 9 * median_flow_obs
        for k in range(len(self.eval_criteria_monthly)):
            eval = self.eval_criteria_monthly[k]
            eval_values_obs = [np.nan] * len(df_obs_monthly)
            eval_values_sim = [np.nan] * len(df_obs_monthly)
            for i, time in enumerate(df_obs_monthly.index):
                m = time.month
                y = time.year
                data_by_month_obs = df_obs[(df_obs.index.month == m) & (df_obs.index.year == y)].copy()
                data_by_month_sim = df_sim[(df_sim.index.month == m) & (df_sim.index.year == y)].copy()
                if not data_by_month_obs.empty:
                    if eval['metric'] == 'Q_mean':
                        eval_values_obs[i] = data_by_month_obs[eval['variable_to_analyze']].mean()
                        eval_values_sim[i] = data_by_month_sim[eval['variable_to_analyze']].mean()
                    if eval['metric'] == 'RR':
                        eval_values_obs[i] = data_by_month_obs[eval['variable_to_analyze']].sum() / \
                                             data_by_month_obs['Rainfall'].sum()
                        eval_values_sim[i] = data_by_month_sim[eval['variable_to_analyze']].sum() / \
                                             data_by_month_obs['Rainfall'].sum()
                    if eval['metric'] == 'high_flow_freq':
                        # https://tosshtoolbox.github.io/TOSSH/matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_frequency.html
                        high_Q_num_obs = sum(data_by_month_obs['Flow'].values > high_flow_obs)
                        eval_values_obs[i] = high_Q_num_obs / len(data_by_month_obs['Flow'])
                        high_Q_num_sim = sum(data_by_month_sim['Flow'].values > high_flow_obs)
                        eval_values_sim[i] = high_Q_num_sim / len(data_by_month_sim['Flow'])
            if k==0:
                data = {eval['metric'] + '_obs': eval_values_obs, eval['metric'] + '_sim': eval_values_sim,
                        'run_id': [nth_run] * len(eval_values_sim)}
                eval_monthly_for_a_run = pd.DataFrame(data, index=df_obs_monthly.index)
            else:
                eval_monthly_for_a_run[eval['metric'] + '_obs'] = eval_values_obs
                eval_monthly_for_a_run[eval['metric'] + '_sim'] = eval_values_sim
            eval_monthly_for_a_run[eval['metric'] + '_bias'] = eval_monthly_for_a_run[eval['metric'] + '_sim'] - eval_monthly_for_a_run[
                eval['metric'] + '_obs']

        return [nth_run, eval_result_for_a_run, eval_monthly_for_a_run]

    def save_results_to_df(self, all_results):

        # ===============================================================
        # Save results from all runs
        # ===============================================================

        ## Store GLUE results (behavioral vs. non-behavioral)
        self.glue_results = [np.nan] * len(all_results)
        self.run_id = [np.nan] * len(all_results)
        for i in range(len(all_results)):
            self.run_id[i] = all_results[i][0]

        ## Store Evaluation metrics for all runs
        self.eval = [np.nan] * len(all_results)
        for i in range(len(all_results)):
            self.eval[i] = all_results[i][1]

        eval_values = np.empty((len(all_results), len(self.eval_names)))
        eval_values[:] = np.nan
        for j in range(len(self.eval)):
            season_idx = -9999
            for i in range(len(self.eval_names)):
                if 'season_transition' in self.eval_names[i]:
                    if season_idx == -9999:
                        season_idx = i
                        season_num = 0
                    eval_values[j][i] = self.eval[j][season_idx][season_num]
                    season_num += 1
                else:
                    eval_values[j][i] = self.eval[j][i]
        self.df_eval = pd.DataFrame(eval_values, index=self.run_id, columns=self.eval_names)

        for i in range(len(all_results)):
            if i==0:
                self.df_eval_mo = all_results[i][2]
            else:
                self.df_eval_mo = pd.concat([self.df_eval_mo, all_results[i][2]])


    def to_csv(self, df_param_to_calibrate=None):
        print('--- Saving data into csv file ---')

        self.df_eval.to_csv(os.path.join(self.out_path, 'post_evaluations.csv'), sep=',', header=True, index=True,
                                 encoding='utf-8', na_rep='nan')
        self.df_eval_mo.to_csv(os.path.join(self.out_path, 'post_evaluations_monthly_metrics.csv'), sep=',', header=True, index=True,
                                 encoding='utf-8', na_rep='nan')



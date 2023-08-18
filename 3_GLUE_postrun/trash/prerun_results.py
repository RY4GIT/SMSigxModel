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

# Global variables
quantiles = [0.10, 0.5, 0.90]

# Global function
def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    # Code from https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy/32216049
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be within [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def triangle_weight(x, a, b, c):
    # a: lowerlim, c:upperlim, b:midpoint
    y = np.full((len(x),), np.nan)
    # for i in range(len(x)):
    y = np.where((x <= a) | (x >= c), 0, y)
    y = np.where((a <= x) & (x <= b), (x - a) / (b - a), y)
    y = np.where((b <= x) & (x <= c), (c - x) / (c - b), y)
    return y


# GLUE object
class GLUE(object):
    def __init__(self, out_path='./', config_path_CFE='', path_GLUE_output='', eval_criteria=dict(), senario_id=9999):

        self.out_path = out_path  # Output folder path
        self.var_names = ["Flow", "Soil Moisture Content"]  # Variables to be analyzed
        self.config_path_CFE = config_path_CFE
        self.path_GLUE_output = path_GLUE_output
        self.eval_criteria = eval_criteria
        self.senario_id = senario_id
        self.out_path_per_senario = os.path.join(out_path, f'senario_{senario_id}')
        
        if not os.path.exists(self.out_path_per_senario):
            os.mkdir(self.out_path_per_senario)
        
    def post_evaluation(self, plot=True):
        # Judge behavioral vs. non-behavioral
        
        # Parameter bounds defined from an excel file
        try:
            df_pri_paras = pd.read_csv(os.path.join(self.path_GLUE_output, 'parameter_priori.csv'), index_col=0)
        except:
            df_pri_paras = pd.read_csv(os.path.join(self.path_GLUE_output, 'parameter_priori.xls'), index_col=0)
        
        try:
            df_eval_metrics = pd.read_csv(os.path.join(self.path_GLUE_output, 'evaluations.csv'), index_col=0)
        except:
            df_eval_metrics = pd.read_csv(os.path.join(self.path_GLUE_output, 'evaluations.xls'), index_col=0)
        
        try:
            df_eval_metrics_monthly = pd.read_csv(os.path.join(self.path_GLUE_output, 'post_evaluations_monthly_metrics.csv'), index_col=0)
        except:
            df_eval_metrics_monthly = pd.read_csv(os.path.join(self.path_GLUE_output, 'post_evaluations_monthly_metrics.xls'), index_col=0)
        
        
        self.eval_names = []
        self.eval_thresh = []
        for i in range(len(self.eval_criteria)):
            if self.eval_criteria[i]["metric"] == 'season_transition':
                self.eval_names.append(
                    f'{self.eval_criteria[i]["metric"]}_dry2wet_s')
                self.eval_names.append(
                    f'{self.eval_criteria[i]["metric"]}_dry2wet_e')
                self.eval_names.append(
                    f'{self.eval_criteria[i]["metric"]}_wet2dry_s')
                self.eval_names.append(
                    f'{self.eval_criteria[i]["metric"]}_wet2dry_e')
                self.eval_thresh.append(self.eval_criteria[i]["threshold"])
                self.eval_thresh.append(self.eval_criteria[i]["threshold"])
                self.eval_thresh.append(self.eval_criteria[i]["threshold"])
                self.eval_thresh.append(self.eval_criteria[i]["threshold"])
            else:
                self.eval_names.append(f'{self.eval_criteria[i]["metric"]} on {self.eval_criteria[i]["variable_to_analyze"]}')
                self.eval_thresh.append(self.eval_criteria[i]["threshold"])

        for i in range(len(self.eval_names)):
            if ('KGE' in self.eval_names[i]) or ('NSE' in self.eval_names[i]):
                df_eval_metrics[self.eval_names[i]+'_Behavioral'] = (df_eval_metrics[self.eval_names[i]] > self.eval_thresh[i])
            elif 'season_transition' in self.eval_names[i]:
                df_eval_metrics[self.eval_names[i]+'_Behavioral'] = (df_eval_metrics[self.eval_names[i]] < self.eval_thresh[i])
        
        glue_judge_column = df_eval_metrics.filter(like='_Behavioral').columns
        df_eval_metrics['GLUE_results_Behavioral'] = df_eval_metrics[glue_judge_column].all(axis='columns')
        
        # Generate random parameters
        self.behavioral_idx = df_eval_metrics[df_eval_metrics['GLUE_results_Behavioral']].index
        
        # Get only behavioral runs 
        self.df_post_paras = df_pri_paras.iloc[self.behavioral_idx].copy()
        self.df_post_eval = df_eval_metrics.iloc[self.behavioral_idx].copy()
        self.df_post_eval_mo = df_eval_metrics_monthly[df_eval_metrics_monthly['run_id'].isin(self.behavioral_idx.values)].copy()
        behavioral_params = len(self.behavioral_idx) * [None]
        for i, run_id in enumerate(self.behavioral_idx):
            behavioral_params[i] = [run_id, self.df_post_paras.loc[run_id]]
            
        self.df_post_eval.to_csv(os.path.join(self.out_path_per_senario, 'post_evaluations.csv'), sep=',', header=True,
                                    index=True, encoding='utf-8', na_rep='nan')
        self.df_post_eval_mo.to_csv(os.path.join(self.out_path_per_senario, 'post_evaluations_monthly_metrics.csv'), sep=',', header=True,
                                    index=True, encoding='utf-8', na_rep='nan')
        
        print('------------- Saved post-evaluation results --------------')
    
        if plot:
            
            nparas = len(df_pri_paras.columns)
            
            # Histogram
            f = plt.figure(figsize=(4 * 4, 4 * 4))
            for i in range(nparas):
                target_para = df_pri_paras.columns[i]
                ax1 = f.add_subplot(4, 4, i + 1)

                df_pri_paras[target_para].plot.hist(bins=10, alpha=0.4, ax=ax1, color="#3182bd", label="Prior")

                if hasattr(self, 'df_post_paras'):
                    self.df_post_paras[target_para].plot.hist(bins=10, alpha=0.8, ax=ax1, color="#3182bd", label="Posterior")

                ax1.set_xlabel(target_para)
                ax1.legend()

                if i != 0:
                    ax1.yaxis.set_visible(False)

            f.savefig(os.path.join(self.out_path_per_senario,'param_dist.png'), dpi=600)
            del f
            print('------------- Saved figures -------------------')

        # Dotty plot
        # Prior vs. posterior parameter distributions
            if hasattr(self, 'df_post_paras'):

                for j in range(len(self.df_post_eval.columns)):
                    f = plt.figure(figsize=(4 * 4, 4 * 4), constrained_layout=True)
                    f.tight_layout()
                    if '_Behavioral' in self.df_post_eval.columns[j]:
                        None
                    else:
                        target_eval = self.df_post_eval.columns[j]
                        for i in range(nparas):
                            target_para = df_pri_paras.columns[i]
                            ax1 = f.add_subplot(4, 4, i + 1)
                            ax1.scatter(self.df_post_paras[target_para], self.df_post_eval[target_eval], alpha=0.5)
                            ax1.tick_params(axis='both', which='major', labelsize=16)
                            ax1.set_xlabel(target_para, fontsize=16)
                            ax1.set_ylabel(target_eval, fontsize=16)

                    # if i !=0:
                    #     ax1.yaxis.set_visible(False)

                    # f.plot()
                    f.savefig(os.path.join(self.out_path_per_senario,"param_dotty_%s.png" % (target_eval)), dpi=600)
                    
                del f
                print('------------- Saved figures -------------------')

        # Dotty plot
        # Parameter interactions for the behavioral runs
            param_interest = ['bb',
                              'satdk',
                              'slop',
                              'satpsi',
                              'smcmax',
                              'wltsmc',
                              'alpha_fc',
                              'lksatfac',
                              'D',
                              'trigger_z_fact',
                              'max_gw_storage',
                              'Cgw',
                              'expon',
                              'refkdt',
                              'K_nash'
                              ]

            if hasattr(self, 'df_post_paras'):
                f = plt.figure(figsize=(20, 20), constrained_layout=True)
                f.tight_layout()
                n_plot = 0
                for i in range(len(param_interest)):
                    for j in range(len(param_interest)):
                        n_plot += 1
                        para0 = param_interest[j]
                        para1 = param_interest[i]
                        ax1 = f.add_subplot(len(param_interest), len(param_interest), n_plot)
                        x = self.df_post_paras[para0]
                        y = self.df_post_paras[para1]
                        ax1.scatter(x, y, alpha=0.5)

                        if i == 0:
                            ax1.xaxis.set_label_position('top')
                            ax1.set_xlabel(para0)
                        if j == 0:
                            ax1.set_ylabel(para1)
                        ax1.tick_params(direction="in")

                f.savefig(os.path.join(self.out_path_per_senario, "param_dotty_interaction.png"), dpi=600)
                del f
                print('------------- Saved figures -------------------')
        
        return behavioral_params


    def reproduce_behavioral_run(self, behavioral_param):

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

            # Merge observed and simulated timeseries
            df = pd.merge_asof(sim, obs, on="Time")
            self.df_timeaxis = df["Time"]

            sim_synced[var_name] = df[var_name + "_x"].copy()
            obs_synced[var_name] = df[var_name + "_y"].copy()

        self.df_behavioral_Q = sim_synced['Flow']
        self.df_behavioral_SM = sim_synced['Soil Moisture Content']

        return [nth_run, self.df_behavioral_Q, self.df_behavioral_SM]


    def calc_uncertainty_bounds(self, all_results, plot=True):
        print('--- Post-processing the be results ---')
        
        # ==================================================
        # One more last run to get synced observed timeseries
        # ==================================================
        myCFE = bmi_cfe.BMI_CFE(self.config_path_CFE)
        myCFE.initialize()
        sim0 = myCFE.run_unit_test(plot=False, warm_up=True)
        obs0 = myCFE.load_unit_test_data()
        obs_synced = pd.DataFrame()

        # Get the results
        for var_name in self.var_names:
            # Get the simulated data
            sim = sim0[["Time", var_name]].copy()
            sim["Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S")  # Works specifically for CFE

            # Get the comparison data
            obs = obs0[["Time", var_name]].copy()
            obs["Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M")  # Works specifically for Mahurangi data

            # Merge observed and simulated timeseries
            df = pd.merge_asof(sim, obs, on="Time")
            self.df_timeaxis = df["Time"]
            obs_synced[var_name] = df[var_name + "_y"].copy()
            
            self.obs_synced = obs_synced
            
        self.df_obs_Q = pd.DataFrame(self.obs_synced["Flow"])
        self.df_obs_Q.set_axis(self.df_timeaxis, axis=0, inplace=True)
        self.df_obs_SM = pd.DataFrame(self.obs_synced["Soil Moisture Content"])
        self.df_obs_SM.set_axis(self.df_timeaxis, axis=0, inplace=True)
            
        # Store Behavioral runs for Flow and soil moisture in dataframe 
        self.run_id_behavioral = []
        for i in range(len(all_results)): 
            self.run_id_behavioral.append(all_results[i][0])
            
        for i in range(len(all_results)):
            if not hasattr(self, 'df_behavioral_Q'):
                self.df_behavioral_Q = all_results[i][1]
                self.df_behavioral_SM = all_results[i][2]
            else:
                self.df_behavioral_Q = pd.concat([self.df_behavioral_Q, all_results[i][1]], axis=1)
                self.df_behavioral_SM = pd.concat([self.df_behavioral_SM, all_results[i][2]], axis=1)

        self.df_behavioral_Q.set_axis(self.run_id_behavioral, axis=1, inplace=True)
        self.df_behavioral_Q.set_axis(self.df_timeaxis, axis=0, inplace=True)
        self.df_behavioral_SM.set_axis(self.run_id_behavioral, axis=1, inplace=True)
        self.df_behavioral_SM.set_axis(self.df_timeaxis, axis=0, inplace=True)

        # Calculate weights
        weight = np.empty((len(self.df_post_eval), len(self.eval_names)))
        j = int(0)
        # Loop for all evaluation metrics
        for i in range(len(self.eval_criteria)):
            if self.eval_criteria[i]['metric'] == "NSE" or self.eval_criteria[i]['metric'] == "KGE":
                # For Nash-Sutcliffe and Kling-Gupta Efficiency scores
                weight[:, j] = ((self.df_post_eval[self.eval_names[j]] - self.eval_criteria[i]['threshold']) / sum(
                    self.df_post_eval[self.eval_names[j]] - self.eval_criteria[i]['threshold'])).to_numpy()
                j += int(1)
            elif self.eval_criteria[i]['metric'] == "season_transition":
                for k in range(4):
                    # For seasonal transition dates
                    weight[:, j + k] = triangle_weight(self.df_post_eval[self.eval_names[j + k]],
                                                       a=-1 * self.eval_criteria[i]['threshold'], b=0,
                                                       c=self.eval_criteria[i]['threshold'])
                j += int(4)
        avg_weight = np.mean(weight, axis=1)

        # Calculate weighted quantile
        for var_name in self.var_names:
            if var_name == "Flow":
                df_behavioral = self.df_behavioral_Q.copy()
                t_len = len(self.df_behavioral_Q)
            elif var_name == "Soil Moisture Content":
                df_behavioral = self.df_behavioral_SM.copy()
                t_len = len(self.df_behavioral_SM)
            np_behavioral = df_behavioral.to_numpy(copy=True)

            # Get weighted quantile
            quantile = np.empty((t_len, len(quantiles)))
            quantile[:] = np.nan
            for t in range(t_len):
                values = np_behavioral[t, :]  # df_behavioral.iloc[[t]].values.flatten()
                quantile[t, :] = weighted_quantile(values=values, quantiles=quantiles, sample_weight=avg_weight,
                                                   values_sorted=False, old_style=False)
            df_simrange = pd.DataFrame(quantile, index=df_behavioral.index, columns=['lowerlim', 'median', 'upperlim'])
            if var_name == "Flow":
                self.df_Q_simrange = df_simrange.copy()
            elif var_name == "Soil Moisture Content":
                self.df_SM_simrange = df_simrange.copy()
        
        if hasattr(self, 'df_Q_simrange'):
            self.df_Q_simrange.to_csv(os.path.join(self.out_path_per_senario,'quantiles_Q.csv'), sep=',', header=True, index=True,
                                      encoding='utf-8', na_rep='nan')
        if hasattr(self, 'df_SM_simrange'):
            self.df_SM_simrange.to_csv(os.path.join(self.out_path_per_senario, 'quantiles_SM.csv'), sep=',', header=True,
                                       index=True, encoding='utf-8', na_rep='nan')

        if plot:
        
            # Plot 
            for var_name in self.var_names:
                if var_name == "Flow":
                    df_simrange = self.df_Q_simrange
                    df_obs = self.df_obs_Q
                    obs_label = 'Observed flow'
                    y_label = 'Total flow [mm/hour]'
                    title = 'Total flow for behavioral runs'
                    fn = 'Q_range.png'
                elif var_name == "Soil Moisture Content":
                    var_name = "Soil Moisture Content"
                    df_simrange = self.df_SM_simrange
                    df_obs = self.df_obs_SM
                    obs_label = 'Observed soil moisture'
                    y_label = 'Volumetric Soil Moisture Content [m^3/m^3]'
                    title = 'Soil moisture for behavioral runs'
                    fn = 'SM_range.png'

                f2 = plt.figure(figsize=(8, 6))
                ax2 = f2.add_subplot()

                df_simrange['lowerlim'].plot(color='black', alpha=0.2, ax=ax2, label=f'{quantiles[0]*100} percentile')
                df_simrange['upperlim'].plot(color='black', alpha=0.2, ax=ax2, label=f'{quantiles[2]*100} percentile')
                df_obs[var_name].plot(color='black', alpha=1, ax=ax2, label=obs_label)
                plt.fill_between(df_simrange.index, df_simrange['upperlim'], df_simrange['lowerlim'],
                                 facecolor='green', alpha=0.2, interpolate=True, label='Predicted range')
                if var_name == "Flow":
                    ax2.set_yscale('log')
                    ax2.set_ylim([1E-07, 1E-01])
                ax2.set_xlabel('Time')
                ax2.set_ylabel(y_label)
                ax2.set_title(title)
                ax2.legend()
                f2.autofmt_xdate()
                f2.savefig(os.path.join(self.out_path_per_senario, fn), dpi=600)



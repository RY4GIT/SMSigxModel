import spotpy
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import sys
sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/libs/cfe/py_cfe")
import cfe
import bmi_cfe
import os
import sys
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import spotpy


sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/libs/SMSig")
from sig_seasontrans import SMSig

# Global variables
quantiles = [0.05, 0.5, 0.95]

# Global function
def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    # Code from https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy/32216049
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
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
    #for i in range(len(x)):
    y = np.where((x<=a) | (x>=c), 0, y)
    y = np.where((a<=x) & (x<=b), (x-a)/(b-a), y)
    y = np.where((b<=x) & (x<=c), (c-x)/(c-b), y)
    return y

# GLUE object
class MyGLUE(object):
    def __init__(self, out_path='./', config_path_CFE='', config_path='', nrun = 1, eval_criteria=dict()):

        self.out_path = out_path # Output folder path
        self.nrun = nrun # Number of runs
        self.var_names = ["Flow", "Soil Moisture Content"] # Variables to be analyzed

        # CFE model instance
        cfe_instance = bmi_cfe.BMI_CFE(config_path_CFE)
        cfe_instance.initialize()
        self.myCFE = cfe_instance

        # Parameter bounds defined from an excel file
        df = pd.read_excel(config_path)
        df_param_to_calibrate = df[df['calibrate'] == 1]
        params = len(df_param_to_calibrate)*[None]
        for i in range(len(params)):
            params[i] = spotpy.parameter.Uniform(
                df_param_to_calibrate['name'][i],
                low=df_param_to_calibrate['lower_bound'][i],
                high=df_param_to_calibrate['upper_bound'][i]
            )
        self.params = params

        # Evaluation criteria (multi-criteria allowed)
        self.eval_criteia = eval_criteria
        print(f"A number of criterion: {len(eval_criteria)}")
        self.eval_names = []
        for i in range(len(eval_criteria)):
            if eval_criteria[i]["metric"]=='season_transition':
                self.eval_names.append(f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (d2w_start)')
                self.eval_names.append(f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (d2w_end)')
                self.eval_names.append(f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (w2d_start)')
                self.eval_names.append(f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (w2d_end)')
            else:
                self.eval_names.append(f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]}')
            print(f'[{i+1}] {eval_criteria[i]["metric"]}-based analysis on {eval_criteria[i]["variable_to_analyze"]}')

        # Initializations
        self.post_rid = []
        self.pri_paras = []
        self.post_paras = []
        self.resulted_totalQ = []
        self.eval = []

    def simulation(self):

        flag_behavioral = 0
        for n in range(self.nrun):

            # ===============================================================
            # Write the randomly-generated parameters to the config json file
            # ===============================================================

            # Generate parameters
            self.sampled = spotpy.parameter.generate(self.params)

            # Get the model config file
            with open(self.myCFE.cfg_file) as data_file:
                self.cfe_cfg = json.load(data_file)

            # Overwrite the model config file
            for i in range(len(self.sampled)):
                if self.sampled[i][1] in ['bb', 'satdk', 'satpsi', 'slop', 'smcmax', 'wltsmc', 'exponent_secondary', 'coeff_secondary', 'D']:
                    self.cfe_cfg["soil_params"][self.sampled[i][1]] = self.sampled[i][0]
                else:
                    self.cfe_cfg[self.sampled[i][1]] = self.sampled[i][0]

            # Save the config file with new parameters
            with open(self.myCFE.cfg_file, 'w') as out_file:
                json.dump(self.cfe_cfg, out_file)

            # Store the prior paramters
            self.pri_paras.append(self.sampled)

            # ===============================================================
            # Actual model run
            # ===============================================================
            self.myCFE.initialize()
            sim0 = self.myCFE.run_unit_test(plot=False)
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
                sim["Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S") # Works specifically for CFE

                # Get the comparison data
                obs = obs0[["Time", var_name]].copy()
                obs["Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M") # Works specifically for Mahurangi data

                # Merge observed and simulated timeseries
                df = pd.merge_asof(sim, obs, on="Time")

                sim_synced[var_name] = df[var_name + "_x"].copy()
                obs_synced[var_name] = df[var_name + "_y"].copy()


            # ===============================================================
            # Evalute the outputs
            # ===============================================================

            # Preparation
            eval_result_for_a_run = []
            behavioral_flag = [False]*len(self.eval_criteia)

            # Loop for all evaluation metrics (multi-criteria).
            # Calculate the metrics and judge behavioral vs. non-behavioral
            for i in range(len(self.eval_criteia)):

                # Nash-Sutcliffe scores
                if self.eval_criteia[i]['metric'] == "NSE":
                    metric_value = spotpy.objectivefunctions.nashsutcliffe(
                        obs_synced[self.eval_criteia[i]['variable_to_analyze']],
                        sim_synced[self.eval_criteia[i]['variable_to_analyze']]
                    )
                    if metric_value > self.eval_criteia[i]['threshold']:
                        behavioral_flag[i] = True

                # Kling-Gupta Efficiency scores
                #　Kling-Gupta efficiencies range from -Inf to 1. Essentially, the closer to 1, the more accurate the model is
                elif self.eval_criteia[i]['metric'] == "KGE":
                    metric_value = spotpy.objectivefunctions.kge(
                        obs_synced[self.eval_criteia[i]['variable_to_analyze']],
                        sim_synced[self.eval_criteia[i]['variable_to_analyze']]
                    )
                    if metric_value > self.eval_criteia[i]['threshold']:
                        behavioral_flag[i] = True

                # Seasonal transition dates
                elif self.eval_criteia[i]['metric'] == "season_transition":
                    # Calculate metrics for OBSERVED timeseries as a baseline performance. Run only once
                    if n == 0:
                        sig_obs = SMSig(
                            ts_time=df["Time"].to_numpy(),
                            ts_value=obs_synced[self.eval_criteia[i]['variable_to_analyze']].to_numpy(),
                            plot_results=False,
                            plot_label="obs"
                        )
                        # sig_obs.detrend() # TODO:debug
                        sig_obs.movmean()
                        t_valley = sig_obs.calc_sinecurve()
                        season_trans_obs = sig_obs.calc_seasontrans(t_valley=t_valley)

                    # Calculate metrics for SIMULATED timeseries
                    sig_sim = SMSig(
                        ts_time=df["Time"].to_numpy(),
                        ts_value=sim_synced[self.eval_criteia[i]['variable_to_analyze']].to_numpy(),
                        plot_results=False,
                        plot_label="sim"
                    )
                    sig_sim.movmean()
                    season_trans_sim = sig_sim.calc_seasontrans(t_valley=t_valley)

                    # Get the deviations in seasonal transition dates between simulated and observed timeseries
                    diff = season_trans_sim - season_trans_obs
                    metric_value = abs(np.nanmean(diff, axis=0))
                    if all(metric_value) < self.eval_criteia[i]['threshold']:
                        behavioral_flag[i] = True

                # Store evaluation metrics for all criteria for one run
                eval_result_for_a_run.append(metric_value)


            # ===============================================================
            # Judge behavioral vs. non-behavioral
            # ===============================================================

            if all(behavioral_flag):
                # If all criteria is TRUE, the model is behavioral
                result_glue = 'Behavioral'

                # Store the behavioral runs
                self.post_rid.append(n) # runid
                self.post_paras.append(self.sampled) # parameters
                self.eval.append(eval_result_for_a_run) # evaluation metrics

                # timeseires
                if flag_behavioral == 0:
                    self.df_behavioral_Q = sim_synced["Flow"]
                    self.df_behavioral_SM = sim_synced["Soil Moisture Content"]
                    flag_behavioral = 1
                else:
                    self.df_behavioral_Q = pd.concat([self.df_behavioral_Q, sim_synced["Flow"]], axis=1)
                    self.df_behavioral_SM = pd.concat([self.df_behavioral_SM, sim_synced["Soil Moisture Content"]], axis=1)

            else:
                # Discard non-behavioral runs
                result_glue = 'Non-behavioral'

            print(f"{n}-th run: {result_glue}")
            print(eval_result_for_a_run)

        # ===============================================================
        # Save results from all runs
        # ===============================================================
        # Get the number of behavioral runs
        n_behavioral = len(self.post_rid)

        # Store PRIOR parameters
        param_names = []
        for i in range(len(self.params)):
            param_names.append(self.pri_paras[0][i][1])
        param_values = np.empty((self.nrun, len(param_names)))
        param_values[:] = np.nan
        for i in range(len(param_names)):
            for j in range(self.nrun):
                param_values[j][i] = self.pri_paras[j][i][0]

        self.df_pri_paras = pd.DataFrame(param_values, columns=param_names)

        # Store POSTERIOR paramters for behavioral runs
        param_names = []
        if n_behavioral!=0:
            # If there is any behavioral parameter
            for i in range(len(self.params)):
                param_names.append(self.post_paras[0][i][1])
            param_values = np.empty((n_behavioral, len(param_names)))
            param_values[:] = np.nan
            for i in range(len(param_names)):
                for j in range(n_behavioral):
                    param_values[j][i] = self.post_paras[j][i][0]
            self.df_post_paras = pd.DataFrame(param_values, index=self.post_rid, columns=param_names)

        # Store Evaluation metrics for behavioral runs
        if n_behavioral != 0:
            eval_values = np.empty((n_behavioral, len(self.eval_names)))
            eval_values[:] = np.nan
            for i in range(len(self.eval_names)):
                for j in range(n_behavioral):
                    if 'season_transition' in self.eval_names[i]:
                        eval_values[j][i] = self.eval[j][0][i]
                    else:
                        eval_values[j][i] = self.eval[j][i]
            self.df_post_eval = pd.DataFrame(eval_values, index=self.post_rid, columns=self.eval_names)

        # Store Simulated timeseries for behavioral runs
        self.df_behavioral_Q.set_axis(self.post_rid, axis=1, inplace=True)
        self.df_behavioral_Q.set_axis(df["Time"], axis=0, inplace=True)
        self.df_behavioral_SM.set_axis(self.post_rid, axis=1, inplace=True)
        self.df_behavioral_SM.set_axis(df["Time"], axis=0, inplace=True)

        # Store Observed timeseries
        self.df_obs_Q = pd.DataFrame(obs_synced["Flow"])
        self.df_obs_Q.set_axis(df["Time"], axis=0, inplace=True)
        self.df_obs_SM = pd.DataFrame(obs_synced["Soil Moisture Content"])
        self.df_obs_SM.set_axis(df["Time"], axis=0, inplace=True)

    def post_process(self):
        # post-process the results
        print('--- Post-processing the simulated results ---')

        # Calculate weights
        weight = np.empty((len(self.df_post_eval), len(self.eval_names)))
        j=int(0)
        # Loop for all evaluation metrics
        for i in range(len(self.eval_criteia)):
            if self.eval_criteia[i]['metric'] == "NSE" or self.eval_criteia[i]['metric'] == "KGE":
                # For Nash-Sutcliffe and Kling-Gupta Efficiency scores
                weight[:,j] = ((self.df_post_eval[self.eval_names[j]] - self.eval_criteia[i]['threshold'])/sum(self.df_post_eval[self.eval_names[j]] - self.eval_criteia[i]['threshold'])).to_numpy()
                j += int(1)
            elif self.eval_criteia[i]['metric'] == "season_transition":
                # For seasonal transition dates
                weight[:, j:j+4] = triangle_weight(self.df_post_eval[self.eval_names[j]], a= -1*self.eval_criteia[i]['threshold'], b=0, c=self.eval_criteia[i]['threshold'])
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
                values = np_behavioral[t,:] # df_behavioral.iloc[[t]].values.flatten()
                quantile[t, :] = weighted_quantile(values=values, quantiles=quantiles, sample_weight=avg_weight,
                                  values_sorted=False, old_style=False)
            df_simrange = pd.DataFrame(quantile, index=df_behavioral.index, columns=['lowerlim', 'median', 'upperlim'])
            if var_name == "Flow":
                self.df_Q_simrange = df_simrange.copy()
            elif var_name == "Soil Moisture Content":
                self.df_SM_simrange = df_simrange.copy()

    def to_csv(self):
        print('--- Saving data into csv file ---')

        # Dump the parameter range to txt file
        file = open(os.path.join(self.out_path, "param_bounds.txt"), "w")
        file.write("%s" % str(self.params))
        file.close

        # Dump the results to csv
        self.df_pri_paras.to_csv(os.path.join(self.out_path, 'paramter_priori.csv'), sep=',', header=True, index=False, encoding='utf-8', na_rep='nan')
        if hasattr(self, 'df_post_paras'):
            self.df_post_paras.to_csv(os.path.join(self.out_path, 'parameter_posterior.csv'), sep=',', header=True, index=False, encoding='utf-8', na_rep='nan')
        if hasattr(self, 'df_post_eval'):
            self.df_post_eval.to_csv(os.path.join(self.out_path, 'evaluations.csv'), sep=',', header=True, index=False, encoding='utf-8', na_rep='nan')
        if hasattr(self, 'df_Q_simrange'):
            self.df_Q_simrange.to_csv(os.path.join(self.out_path, 'quantiles_Q.csv'), sep=',', header=True, index=False, encoding='utf-8', na_rep='nan')
        if hasattr(self, 'df_SM_simrange'):
            self.df_SM_simrange.to_csv(os.path.join(self.out_path, 'quantiles_SM.csv'), sep=',', header=True, index=False, encoding='utf-8', na_rep='nan')

    def plot(self, plot_type=None):
        # Plot the results
        print('--- Saving the plots ---')

        # Histogram
        # Prior vs. posterior parameter distributions
        if plot_type == "param_hist":
            nparas = len(self.df_pri_paras.columns)
            f = plt.figure(figsize=(4*4, 4*4))

            for i in range(nparas):
                target_para = self.df_pri_paras.columns[i]
                ax1 = f.add_subplot(4, 4, i+1)

                self.df_pri_paras[target_para].plot.hist(bins=10, alpha=0.4, ax=ax1, color="#3182bd", label="Prior")

                if hasattr(self, 'df_post_paras'):
                    self.df_post_paras[target_para].plot.hist(bins=10, alpha=0.8, ax=ax1, color="#3182bd", label="Posterior")

                ax1.set_xlabel(target_para)
                ax1.legend()

                if i !=0:
                    ax1.yaxis.set_visible(False)

            # f.plot()

            f.savefig(os.path.join(self.out_path, 'param_dist.png'), dpi=600)


        # Dotty plot
        # Prior vs. posterior parameter distributions
        if plot_type == "dotty":
            if hasattr(self, 'df_post_paras'):

                nparas = len(self.df_pri_paras.columns)

                for j in range(len(self.df_post_eval.columns)):
                    f = plt.figure(figsize=(4 * 4, 4 * 4), constrained_layout=True)
                    f.tight_layout()
                    target_eval = self.df_post_eval.columns[j]
                    for i in range(nparas):
                        target_para = self.df_pri_paras.columns[i]
                        ax1 = f.add_subplot(4, 4, i+1)
                        plt.scatter(self.df_post_paras[target_para], self.df_post_eval[target_eval])

                        ax1.set_xlabel(target_para)
                        ax1.set_ylabel(target_eval)

                    # if i !=0:
                    #     ax1.yaxis.set_visible(False)

                    # f.plot()
                    f.savefig(os.path.join(self.out_path, "param_dotty_%s.png" % (target_eval)), dpi=600)

        # Dotty plot
        # Parameter interactions for the behavioral runs
        if plot_type == "dotty_interaction":

            param_interset = ['bb',
                      'satdk',
                      'satpsi',
                      'slop',
                      'smcmax',
                      'wltsmc',
                      'refkdt',
                      'trigger_z_m',
                      'fc_atm_press_fraction'
                      ]

            if hasattr(self, 'df_post_paras'):
                    f = plt.figure(figsize=(16, 16), constrained_layout=True)
                    f.tight_layout()
                    n_plot = 0
                    for i in range(len(param_interset)):
                        for j in range(len(param_interset)):
                            n_plot += 1
                            para0 = param_interset[i]
                            para1 = param_interset[j]
                            ax1 = f.add_subplot(len(param_interset), len(param_interset), n_plot)
                            x = self.df_post_paras[para0]
                            y = self.df_post_paras[para1]
                            plt.scatter(x, y)

                            ax1.set_xlabel(para0)
                            ax1.set_ylabel(para1)
                            ax1.tick_params(direction="in")

                            # if i !=0:
                            #     ax1.yaxis.set_visible(False)

                    # f.plot()
                    f.savefig(os.path.join(self.out_path, "param_dotty_interaction.png"), dpi=600)

        # Time series of data
        # Flow and Soil Moisture Content
        if plot_type == "timeseries":
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

                df_simrange['lowerlim'].plot(color='black', alpha=0.2, ax=ax2, label='_Hidden')
                df_simrange['upperlim'].plot(color='black', alpha=0.2, ax=ax2, label='_Hidden')
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
                f2.savefig(os.path.join(self.out_path, fn), dpi=600)


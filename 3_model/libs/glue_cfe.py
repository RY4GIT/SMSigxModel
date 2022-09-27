import spotpy
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import cfe
import os
import sys
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

# Global variables
KGE_Q_thresh = -0.41 # threshold value
KGE_SM_thresh = 0.5 # threshold value
NSE_Q_thresh = -1000
seasontrans_thresh = 35 # seasonal transition threshold
quantiles = [0.05, 0.5, 0.95] # quantiles

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/SMSig")
from sig_seasontrans import SMSig

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
    def __init__(self, cfe_input, out_path='./', nrun = 1, calib_case=1):

        self.out_path = out_path
        self.myCFE = cfe_input

        # define parameter bounds
        self.params = [
            spotpy.parameter.Uniform('bb', low=2, high=15),
            spotpy.parameter.Uniform('satdk', low=0, high=1),
            spotpy.parameter.Uniform('slop', low=0, high=1),
            spotpy.parameter.Uniform('satpsi', low=0.02, high=0.78),
            spotpy.parameter.Uniform('smcmax', low=0.33, high=1.0),
            spotpy.parameter.Uniform('wltsmc', low=0.0, high=0.57),
            spotpy.parameter.Uniform('exponent_secondary', low=0, high=8),
            spotpy.parameter.Uniform('max_gw_storage', low=0.1, high=250),
            spotpy.parameter.Uniform('Cgw', low=0.01, high=10),
            spotpy.parameter.Uniform('expon', low=0, high=8),
            spotpy.parameter.Uniform('K_lf', low=0, high=50),
            spotpy.parameter.Uniform('K_nash', low=0, high=30),
            spotpy.parameter.Uniform('trigger_z_m', low=0.01, high=0.87),
            spotpy.parameter.Uniform('fc_atm_press_fraction', low=0.10, high=0.33),
            spotpy.parameter.Uniform('refkdt', low=0.1, high=4)
        ]

        # define the number of iterations
        self.nrun = nrun

        # initialization
        self.post_rid = []
        self.pri_paras = []
        self.post_paras = []
        self.resulted_totalQ = []
        self.eval = []

        # calibration cases
        self.calib_case = calib_case
        self.var_names = ["Flow", "Soil Moisture Content"]

        self.calib_case = calib_case
        if calib_case == 1:
            print('NSE-based calib on Q')
            self.eval_names = ['NSE_Q']
        elif calib_case == 2:
            print('KGE-based calib on Q')
            self.eval_names = ['KGE_Q']
        elif calib_case == 3:
            print('KGE-based calib on SM')
            self.eval_names = ['KGE_SM']
        elif calib_case == 4:
            print('KGE-based calib on Q&SM')
            self.eval_names = ['KGE_Q', 'KGE_SM']
        elif calib_case == 5:
            print('Seasonsig-based calib on SM')
            self.eval_names = ['d2w_start', 'd2w_end', 'w2d_start', 'w2d_end']
        elif calib_case == 6:
            print('KGE+Seasonsig calib on Q&SM')
            self.eval_names = ['KGE_Q', 'KGE_SM', 'd2w_start', 'd2w_end', 'w2d_start', 'w2d_end']


    def simulation(self):
    # run the simulation
        # print('--- Running CFE model ---')
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

            # ===============================================================
            # Actual model run
            # ===============================================================
            self.myCFE.initialize()
            sim0 = self.myCFE.run_unit_test(plot=False)
            obs0 = self.myCFE.load_unit_test_data()

            # ===============================================================
            # Get the simulated and observed data & evaluate
            # ===============================================================
            for var_name in self.var_names:
                # Get the simulated data
                sim = sim0[["Time", var_name]]
                sim["Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S") # Works specifically for CFE

                # Get the comparison data
                obs = obs0[["Time", var_name]]
                obs["Time"] = pd.to_datetime(obs["Time"], format="%d-%b-%Y %H:%M:%S") # Works specifically for Mahurangi data

                # Merge observed and simulated timeseries
                df = pd.merge_asof(sim, obs, on="Time")

                sim_synced = df[var_name + "_x"]
                obs_synced = df[var_name + "_y"]

                # Model evaluators
                # KGE for both streamflow and SM data
                if "KGE_Q" in self.eval_names or "KGE_SM" in self.eval_names:
                    KGE = spotpy.objectivefunctions.kge(obs_synced, sim_synced) #kge_non_parametric(obs_synced, sim_synced)
                if "NSE_Q" in self.eval_names or "NSE_SM" in self.eval_names:
                    NSE = spotpy.objectivefunctions.nashsutcliffe(obs_synced, sim_synced)

                # Seasonsig for SM data
                if var_name == "Soil Moisture Content" and 'd2w_start' in self.eval_names:
                    # Observed. Run only once
                    if n == 0:
                        sig_obs = SMSig(ts_time=df["Time"].to_numpy(), ts_value=obs_synced.to_numpy(), plot_results=False, plot_label="obs")
                        # sig_obs.detrend() # TODO:debug
                        sig_obs.movmean()
                        t_valley = sig_obs.calc_sinecurve()
                        season_trans_obs = sig_obs.calc_seasontrans(t_valley=t_valley)

                    # simulated
                    sig_sim = SMSig(ts_time=df["Time"].to_numpy(), ts_value=sim_synced.to_numpy(), plot_results=False, plot_label="sim")
                    sig_sim.movmean()
                    season_trans_sim = sig_sim.calc_seasontrans(t_valley=t_valley)

                    diff = season_trans_sim - season_trans_obs
                    diff_avg = abs(np.nanmean(diff, axis=0))

                # Get the variables
                if var_name == "Flow":
                    if "KGE_Q" in self.eval_names:
                        KGE_Q = KGE
                    if "NSE_Q" in self.eval_names:
                        NSE_Q = NSE
                    sim_Q_synced = sim_synced
                    obs_Q_synced = obs_synced
                elif var_name == "Soil Moisture Content":
                    if "KGE_SM" in self.eval_names:
                        KGE_SM = KGE
                    if "NSE_SM" in self.eval_names:
                        NSE_SM = NSE
                    sim_SM_synced = sim_synced
                    obs_SM_synced = obs_synced

                del sim, obs

            del sim0, obs0

            # ===============================================================
            # Judge behavioral vs. non-behavioral
            # ===============================================================
            # Store all the prior parameters
            self.pri_paras.append(self.sampled)
            if self.calib_case == 1:
                behavioral_condition = NSE_Q > NSE_Q_thresh
                eval_array = [NSE_Q]
            elif self.calib_case == 2:
                behavioral_condition = KGE_Q > KGE_Q_thresh
                eval_array = [KGE_Q]
            elif self.calib_case == 3:
                behavioral_condition = KGE_SM > KGE_SM_thresh
                eval_array = [KGE_SM]
            elif self.calib_case == 4:
                behavioral_condition = KGE_Q > KGE_Q_thresh and KGE_SM > KGE_SM_thresh
                eval_array = [KGE_Q, KGE_SM]
            elif self.calib_case == 5:
                behavioral_condition = all(diff_avg < seasontrans_thresh)
                eval_array = [diff_avg[0], diff_avg[1], diff_avg[2], diff_avg[3]]
            elif self.calib_case == 6:
                behavioral_condition = KGE_Q > KGE_Q_thresh and KGE_SM > KGE_SM_thresh and all(diff_avg < seasontrans_thresh)
                eval_array = [KGE_Q, KGE_SM, diff_avg[0], diff_avg[1], diff_avg[2], diff_avg[3]]


            if behavioral_condition:
                print_results = 'Behavioral'
                # This is the threshold conditions
                    # KGE must be above thresholds
                    # Average transition dates are less than threshold
                # Store the behavioral runs
                self.post_rid.append(n) #runid
                self.post_paras.append(self.sampled) #parameters
                if flag_behavioral == 0:
                    if "Flow" in self.var_names:
                        self.df_Q_behavioral = sim_Q_synced #timeseires
                    if "Soil Moisture Content" in self.var_names:
                        self.df_SM_behavioral = sim_SM_synced
                else:
                    if "Flow" in self.var_names:
                        self.df_Q_behavioral = pd.concat([self.df_Q_behavioral, sim_Q_synced], axis=1)
                    if "Soil Moisture Content" in self.var_names:
                        self.df_SM_behavioral = pd.concat([self.df_SM_behavioral, sim_SM_synced], axis=1)
                flag_behavioral = 1
                self.eval.append(eval_array)
            else:
                # Discard non-behavioral runs
                print_results = 'Non-behavioral'

            print('Case{}: {}-th run'.format(self.calib_case, n))
            print(print_results)
            print(eval_array)

        # ===============================================================
        # Save results in Dataframe
        # ===============================================================
        # Get the number of behavioral runs
        n_behavioral = len(self.post_rid)

        # Posterior paramters
        param_names = []
        if len(self.post_paras) !=0:
            for i in range(len(self.params)):
                param_names.append(self.post_paras[0][i][1])

            param_values = np.empty((n_behavioral, len(param_names)))
            param_values[:] = np.nan
            for i in range(len(param_names)):
                for j in range(n_behavioral):
                    param_values[j][i] = self.post_paras[j][i][0]
            self.df_post_paras = pd.DataFrame(param_values, index=self.post_rid, columns=param_names)

        # Prior parameters
        param_names = []
        for i in range(len(self.params)):
            param_names.append(self.pri_paras[0][i][1])

        param_values = np.empty((self.nrun, len(param_names)))
        param_values[:] = np.nan
        for i in range(len(param_names)):
            for j in range(self.nrun):
                param_values[j][i] = self.pri_paras[j][i][0]

        self.df_pri_paras = pd.DataFrame(param_values, columns=param_names)

        # Evaluation metrics
        if len(self.eval) != 0:
            eval_values = np.empty((n_behavioral, len(self.eval_names)))
            eval_values[:] = np.nan
            for i in range(len(self.eval_names)):
                for j in range(n_behavioral):
                    eval_values[j][i] = self.eval[j][i]
            self.df_post_eval = pd.DataFrame(eval_values, index=self.post_rid, columns=self.eval_names)

        # Simulated
        # Total flow
        if hasattr(self, 'df_Q_behavioral'):
            self.df_Q_behavioral.set_axis(self.post_rid, axis=1, inplace=True)
            self.df_Q_behavioral.set_axis(df["Time"], axis=0, inplace=True)
        # Soil moisture
        if hasattr(self, 'df_SM_behavioral'):
            self.df_SM_behavioral.set_axis(self.post_rid, axis=1, inplace=True)
            self.df_SM_behavioral.set_axis(df["Time"], axis=0, inplace=True)

        # Observed
        # Total flow
        if "Flow" in self.var_names:
            self.df_Q_obs = pd.DataFrame(obs_Q_synced)
            self.df_Q_obs.set_axis(df["Time"], axis=0, inplace=True)

        # Soil moisture
        if "Soil Moisture Content" in self.var_names:
            self.df_SM_obs = pd.DataFrame(obs_SM_synced)
            self.df_SM_obs.set_axis(df["Time"], axis=0, inplace=True)

    def post_process(self):
        # post-process the results
        print('--- Post-processing the simulated results ---')

        if hasattr(self, 'df_Q_behavioral') or hasattr(self, 'df_SM_behavioral'):

            # Get empirical weight
            if "KGE_Q" in self.eval_names:
                KGE_Q_weight = (self.df_post_eval["KGE_Q"] - KGE_Q_thresh)/sum( (self.df_post_eval["KGE_Q"] - KGE_Q_thresh))
            if "KGE_SM" in self.eval_names:
                KGE_SM_weight = (self.df_post_eval["KGE_SM"] - KGE_SM_thresh) / sum((self.df_post_eval["KGE_SM"] - KGE_SM_thresh))
            if "NSE_Q" in self.eval_names:
                NSE_Q_weight = (self.df_post_eval["NSE_Q"] - NSE_Q_thresh) / sum((self.df_post_eval["NSE_Q"] - NSE_Q_thresh))
            if "d2w_start" in self.eval_names:
                Seasonsig_weight_d2ws = triangle_weight(self.df_post_eval['d2w_start'], a= -1*seasontrans_thresh, b=0, c=seasontrans_thresh)
                Seasonsig_weight_d2we = triangle_weight(self.df_post_eval['d2w_end'], a=-1 * seasontrans_thresh, b=0, c=seasontrans_thresh)
                Seasonsig_weight_w2ds = triangle_weight(self.df_post_eval['w2d_start'], a=-1 * seasontrans_thresh, b=0, c=seasontrans_thresh)
                Seasonsig_weight_w2de = triangle_weight(self.df_post_eval['w2d_end'], a=-1 * seasontrans_thresh, b=0, c=seasontrans_thresh)

            # Get composite weight
            if self.calib_case == 1:
                weight = NSE_Q_weight
            elif self.calib_case == 2:
                weight = KGE_Q_weight
            elif self.calib_case == 3:
                weight = KGE_SM_weight
            elif self.calib_case == 4:
                weight = (KGE_Q_weight.values + KGE_SM_weight.values )/2
            elif self.calib_case == 5:
                weight = (Seasonsig_weight_d2ws + Seasonsig_weight_d2we + Seasonsig_weight_w2ds + Seasonsig_weight_w2de) / 4
            elif self.calib_case == 6:
                weight = (KGE_Q_weight.values + KGE_SM_weight.values + Seasonsig_weight_d2ws + Seasonsig_weight_d2we + Seasonsig_weight_w2ds + Seasonsig_weight_w2de) / 6

            # Calculate weighted quantile
            for var_name in self.var_names:
                if var_name == "Flow":
                    df_behavioral = self.df_Q_behavioral
                    t_len = len(self.df_Q_behavioral)
                elif var_name == "Soil Moisture Content":
                    df_behavioral = self.df_SM_behavioral
                    t_len = len(self.df_SM_behavioral)
                np_behavioral = df_behavioral.to_numpy(copy=True)

                # Get weighted quantile
                quantile = np.empty((t_len, len(quantiles)))
                quantile[:] = np.nan
                for t in range(t_len):
                    values = np_behavioral[t,:] # df_behavioral.iloc[[t]].values.flatten()
                    quantile[t, :] = weighted_quantile(values=values, quantiles=quantiles, sample_weight=weight,
                                      values_sorted=False, old_style=False)
                df_simrange = pd.DataFrame(quantile, index=df_behavioral.index, columns=['lowerlim', 'median', 'upperlim'])
                if var_name == "Flow":
                    self.df_Q_simrange = df_simrange
                elif var_name == "Soil Moisture Content":
                    self.df_SM_simrange = df_simrange

    def to_csv(self):
        print('--- Saving data into csv file ---')

        # dump the parameter range to txt file
        file = open(os.path.join(self.out_path, "param_bounds.txt"), "w")
        file.write("%s" % str(self.params))
        file.close

        # save the results to csv
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
                    f.savefig(os.path.join(self.out_path, "param_dotty_interaction_case%s.png" % (self.calib_case)), dpi=600)

        # Time series of data
        if plot_type == "timeseries":
            if hasattr(self, 'df_Q_simrange') or hasattr(self, 'df_SM_simrange'):
                for var_name in self.var_names:
                    if var_name == "Flow":
                        df_simrange = self.df_Q_simrange
                        df_obs = self.df_Q_obs
                        obs_label = 'Observed flow'
                        y_label = 'Total flow [mm/hour]'
                        title = 'Total flow for behavioral runs'
                        fn = 'Q_range.png'
                    elif var_name == "Soil Moisture Content":
                        var_name = "Soil Moisture Content"
                        df_simrange = self.df_SM_simrange
                        df_obs = self.df_SM_obs
                        obs_label = 'Observed soil moisture'
                        y_label = 'Volmetric Soil Moisture Content [m^3/m^3]'
                        title = 'Soil moisture for behavioral runs'
                        fn = 'SM_range.png'

                    f2 = plt.figure(figsize=(8, 6))
                    ax2 = f2.add_subplot()

                    df_simrange['lowerlim'].plot(color='black', alpha=0.2, ax=ax2, label='_Hidden')
                    df_simrange['upperlim'].plot(color='black', alpha=0.2, ax=ax2, label='_Hidden')
                    df_obs[var_name + "_y"].plot(color='black', alpha=1, ax=ax2, label=obs_label)
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


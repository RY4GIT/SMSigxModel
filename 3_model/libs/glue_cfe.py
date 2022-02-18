import spotpy
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import cfe
import os
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

# Global variables
var_names = ["Total Discharge", "Soil Moisture Content"] # variables of interst
eval_names = ['KGE_Q', 'KGE_SM'] # evaluators
KGE_Q_thresh = -10 # threshold value
KGE_SM_thresh = -10 # threshold value
quantiles = [0.05, 0.5, 0.95] # quantiles

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
    y = np.where((b<=x) & (x<=c), (b-x)/(c-b), y)
    return y

# GLUE object
class MyGLUE(object):
    def __init__(self, cfe_input, out_path='./', obj_func=None):

        self.obj_func = False
        self.out_path = out_path
        self.myCFE = cfe_input

        # define parameter bounds
        self.params = [
            spotpy.parameter.Uniform('smcmax', low=0.1, high=0.7),
            spotpy.parameter.Uniform('wltsmc', low=0.0, high=0.6)
        ]

        # define the number of iterations
        self.nrun = 3

        # initialization
        self.post_rid = []
        self.pri_paras = []
        self.post_paras = []
        self.resulted_totalQ = []
        self.eval = []

    def simulation(self):
    # run the simulation
        print('--- Running CFE model ---')
        for n in range(self.nrun):

            print('{}-th run'.format(n))

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
                if self.sampled[i][1] in ['depth', 'bb', 'mult', 'satdk', 'satpsi', 'slop', 'smcmax', 'wltsmc', 'D']:
                    self.cfe_cfg["soil_params"][self.sampled[i][1]] = self.sampled[i][0]
                else:
                    self.cfe_cfg[self.sampled[i][1]] = self.sampled[i][0]

            # Save the config file with new parameters
            with open(self.myCFE.cfg_file, 'w') as out_file:
                json.dump(self.cfe_cfg, out_file)


            # ===============================================================
            # Actual model run
            # ===============================================================
            sim0 = self.myCFE.run_unit_test(plot=False)
            obs0 = self.myCFE.load_unit_test_data()

            # ===============================================================
            # Get the simulated and observed data & evaluate
            # ===============================================================
            for var_name in var_names:
                # Get the simulated data
                sim = sim0[["Time", var_name]]
                sim["Time"] = pd.to_datetime(sim["Time"])

                # Get the comparison data
                obs = obs0[["Time", var_name]]
                obs["Time"] = pd.to_datetime(obs["Time"])

                # Merge observed and simulated timeseries
                df = pd.merge_asof(sim, obs, on="Time")
                sim_synced = df[var_name + "_x"]
                obs_synced = df[var_name + "_y"]

                # Model evaluators
                KGE = spotpy.objectivefunctions.kge_non_parametric(obs_synced, sim_synced)

                if var_name == "Total Discharge":
                    KGE_Q = KGE
                    sim_Q_synced = sim_synced
                    obs_Q_synced = obs_synced
                elif var_name == "Soil Moisture Content":
                    KGE_SM = KGE
                    sim_SM_synced = sim_synced
                    obs_SM_synced = obs_synced

            # ===============================================================
            # Judge behavioral vs. non-behavioral
            # ===============================================================
            # Store all the prior parameters
            self.pri_paras.append(self.sampled)
            if KGE_Q > KGE_Q_thresh and KGE_SM > KGE_SM_thresh: # This is the threshold conditions
                # Store the behavioral runs
                self.post_rid.append(n) #runid
                self.post_paras.append(self.sampled) #parameters
                if n == 0:
                    self.df_Q_behavioral = sim_Q_synced #timeseires
                    self.df_SM_behavioral = sim_SM_synced
                else:
                    self.df_Q_behavioral = pd.concat([self.df_Q_behavioral, sim_Q_synced], axis=1)
                    self.df_SM_behavioral = pd.concat([self.df_SM_behavioral, sim_SM_synced], axis=1)
                self.eval.append([KGE_Q, KGE_SM])
            else:
                # Discard non-behavioral runs
                None

        # ===============================================================
        # Save results in Dataframe
        # ===============================================================
        # Get the number of behavioral runs
        n_behavioral = len(self.post_rid)

        # Posterior paramters
        param_names = []
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
        eval_values = np.empty((n_behavioral, len(eval_names)))
        eval_values[:] = np.nan
        for i in range(len(eval_names)):
            for j in range(n_behavioral):
                eval_values[j][i] = self.eval[j][i]
        self.df_post_eval = pd.DataFrame(eval_values, index=self.post_rid, columns=eval_names)

        # Simulated
        # Total flow
        self.df_Q_behavioral.set_axis(self.post_rid, axis=1, inplace=True)
        self.df_Q_behavioral.set_axis(df["Time"], axis=0, inplace=True)
        # Soil moisture
        self.df_SM_behavioral.set_axis(self.post_rid, axis=1, inplace=True)
        self.df_SM_behavioral.set_axis(df["Time"], axis=0, inplace=True)

        # Observed
        # Total flow
        self.df_Q_obs = pd.DataFrame(obs_Q_synced)
        self.df_Q_obs.set_axis(df["Time"], axis=0, inplace=True)

        # Soil moisture
        self.df_SM_obs = pd.DataFrame(obs_SM_synced)
        self.df_SM_obs.set_axis(df["Time"], axis=0, inplace=True)
        # TODO: Check observed values

    def post_process(self):
        # post-process the results
        print('--- Post-processing the simulated results ---')

        # Initialize
        t_len = len(self.df_Q_behavioral)

        # Get empirical weight
        KGE_Q_weight = (self.df_post_eval["KGE_Q"] - KGE_Q_thresh)/sum( (self.df_post_eval["KGE_Q"] - KGE_Q_thresh))
        KGE_SM_weight = (self.df_post_eval["KGE_SM"] - KGE_SM_thresh) / sum((self.df_post_eval["KGE_SM"] - KGE_SM_thresh))
        # Get triangular weight
        # TODO: Get triangular weight

        # Get composite weight
        weight = KGE_Q_weight.values + KGE_SM_weight.values

        # Calculate weighted quantile
        for var_name in var_names:
            if var_name == "Total Discharge":
                df_behavioral = self.df_Q_behavioral
            elif var_name == "Soil Moisture Content":
                df_behavioral = self.df_SM_behavioral
            # Get weighted quantile
            quantile = np.empty((t_len, len(quantiles)))
            quantile[:] = np.nan
            for t in range(t_len):
                values = df_behavioral.iloc[[t]].values.flatten()
                quantile[t, :] = weighted_quantile(values=values, quantiles=quantiles, sample_weight=weight,
                                  values_sorted=False, old_style=False)
            df_simrange = pd.DataFrame(quantile, index=self.df_Q_behavioral.index, columns=['lowerlim', 'median', 'upperlim'])
            if var_name == "Total Discharge":
                self.df_Q_simrange = df_simrange
            elif var_name == "Soil Moisture Content":
                self.df_SM_simrange = df_simrange

    def to_csv(self):
        print('--- Saving data into csv file ---')

        # save the results to csv
        self.df_post_paras.to_csv(os.path.join(self.out_path, 'parameter_posterior.csv'))
        self.df_pri_paras.to_csv(os.path.join(self.out_path, 'paramter_priori.csv'))
        # self.df_Q_behavioral.to_csv(os.path.join(self.out_path, 'posterior_ts_Q.csv'))
        # self.df_SM_behavioral.to_csv(os.path.join(self.out_path, 'posterior_ts_SM.csv'))
        self.df_post_eval.to_csv(os.path.join(self.out_path, 'evaluations.csv'))
        self.df_Q_simrange.to_csv(os.path.join(self.out_path, 'quantiles_Q.csv'))
        self.df_SM_simrange.to_csv(os.path.join(self.out_path, 'quantiles_SM.csv'))

    def plot(self):
        # Plot the results
        print('--- Saving the plots ---')

        # Prior vs. posterior parameter distributions
        nparas = len(self.df_pri_paras.columns)
        f = plt.figure(figsize=(4*nparas, 4))
        for i in range(nparas):
            target_para = self.df_pri_paras.columns[i]
            ax1 = f.add_subplot(1, nparas, i+1)
            self.df_pri_paras[target_para].plot.hist(bins=10, alpha=0.4, ax=ax1, color="#3182bd", label="Prior")
            self.df_post_paras[target_para].plot.hist(bins=10, alpha=0.8, ax=ax1, color="#3182bd", label="Posterior")
            ax1.set_xlabel(target_para)
            ax1.legend()
            if i !=0:
                ax1.yaxis.set_visible(False)
        f.savefig(os.path.join(self.out_path, 'param_dist.png'), dpi=600)

        # Time series of data
        for var_name in var_names:
            if var_name == "Total Discharge":
                df_simrange = self.df_Q_simrange
                df_obs = self.df_Q_obs
                obs_label = 'Observed discharge'
                y_label = 'Total flow [mm/hour]'
                title = 'Total flow for behavioral runs'
                fn = 'Q_range.png'
            elif var_name == "Soil Moisture Content":
                df_simrange = self.df_SM_simrange
                df_obs = self.df_SM_obs
                obs_label = 'Observed soil moisture'
                y_label = 'Volmetric Soil Moisture Content [m^3/m^3]'
                title = 'Soil moisture for behavioral runs'
                fn = 'SM_range.png'

            f2 = plt.figure(figsize=(8,6))
            ax2 = f2.add_subplot()
            df_simrange['lowerlim'].plot(color='black',alpha=0.2, ax= ax2,  label='_Hidden')
            df_simrange['upperlim'].plot(color='black', alpha=0.2, ax=ax2,  label='_Hidden')
            df_obs.plot(color='black', alpha=1, ax=ax2, label=obs_label)
            plt.fill_between(df_simrange.index, df_simrange['upperlim'], df_simrange['lowerlim'],
                             facecolor='green', alpha=0.2, interpolate=True, label='Predicted range')
            ax2.set_xlabel('Time')
            ax2.set_ylabel(y_label)
            ax2.set_title(title)
            ax2.legend()
            f2.savefig(os.path.join(self.out_path, fn), dpi=600)



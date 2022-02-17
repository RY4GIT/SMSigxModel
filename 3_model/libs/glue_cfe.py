import spotpy
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import cfe
import os
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

out_file_path = '../4_out/Mahurangi/'

class MyGLUE(object):
    def __init__(self, cfe_input, obj_func=None):

        self.obj_func = False

        # Model
        self.myCFE = cfe_input

        # define parameter bounds
        self.params = [
            spotpy.parameter.Uniform('smcmax', low=0.1, high=0.7),
            spotpy.parameter.Uniform('wltsmc', low=0.0, high=0.6)
        ]

        self.nrun = 3

        # initialization
        self.post_rid = []
        self.pri_paras = []
        self.post_paras = []
        self.resulted_totalflow = []
        self.eval = []

    def simulation(self):

        for n in range(self.nrun):

            print('{}-th run'.format(n))

            # ===============================================================
            # Write the randomly-generated parameters to the config json file
            # ===============================================================

            self.sampled = spotpy.parameter.generate(self.params)

            with open(self.myCFE.cfg_file) as data_file:
                self.cfe_cfg = json.load(data_file)

            for i in range(len(self.sampled)):
                if self.sampled[i][1] in ['depth', 'bb', 'mult', 'satdk', 'satpsi', 'slop', 'smcmax', 'wltsmc', 'D']:
                    self.cfe_cfg["soil_params"][self.sampled[i][1]] = self.sampled[i][0]
                else:
                    self.cfe_cfg[self.sampled[i][1]] = self.sampled[i][0]

            with open(self.myCFE.cfg_file, 'w') as out_file:
                json.dump(self.cfe_cfg, out_file)

            # ===============================================================
            # Actual model run & get the simulated discharge
            # ===============================================================
            data = self.myCFE.run_unit_test(plot=False)
            sim = data[["Time", "Total Discharge"]]
            sim["Time"] = pd.to_datetime(sim["Time"])
            # return sim

            # ===============================================================
            # Get the observed discharge & evaluate
            # ===============================================================
            # Get the comparison data
            data = self.myCFE.load_unit_test_data()
            obs = data[["Time", "Total Discharge"]]
            obs["Time"] = pd.to_datetime(obs["Time"])
            # return obs

            # Merge observed and simulated timeseries
            df = pd.merge_asof(sim, obs, on = "Time")
            sim_synced = df["Total Discharge_x"]
            obs_synced = df["Total Discharge_y"]

            # Model evaluators
            eval_names = ['KGE']
            like1 = spotpy.objectivefunctions.kge_non_parametric(obs_synced, sim_synced)
            # TODO: add soil moisture, as well as seasonal transition

            # Choose a behavioral threshold and store the good parameters
            self.pri_paras.append(self.sampled)
            if like1 > 0: # This is the threshold conditions
                self.post_rid.append(n)
                self.post_paras.append(self.sampled)
                if n == 0:
                    self.resulted_totalflow = sim_synced
                else:
                    self.resulted_totalflow = pd.concat([self.resulted_totalflow, sim_synced], axis=1)
                self.eval.append([like1])
            else:
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

        # Total flow
        self.df_post_flow = self.resulted_totalflow
        self.df_post_flow.set_axis(self.post_rid, axis=1, inplace=True)
        self.df_post_flow.set_axis(df["Time"], axis=0, inplace=True)

        # Evaluation metrics
        eval_values = np.empty((n_behavioral, len(eval_names)))
        eval_values[:] = np.nan
        for i in range(len(eval_names)):
            for j in range(n_behavioral):
                eval_values[j][i] = self.eval[j][i]
        self.df_post_eval = pd.DataFrame(eval_values, index=self.post_rid, columns=eval_names)

    def post_process(self):
        # Get weighted quantile
        None


    def to_csv(self):
        self.df_post_paras.to_csv(os.path.join(out_file_path, 'posterio_parameter.csv'))
        self.df_pri_paras.to_csv(os.path.join(out_file_path, 'priori_parameter.csv'))
        self.df_post_flow.to_csv(os.path.join(out_file_path, 'posterior_flows.csv'))
        self.df_post_eval.to_csv(os.path.join(out_file_path, 'evaluations.csv'))

    def plot(self):
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
        f.savefig(os.path.join(out_file_path, 'param_dist.png'), dpi=600)


        # Total flow
        # ax2 = plt.figure(figsize=(10,10))
        f2 = plt.figure(figsize=(8,6))
        ax2 = f2.add_subplot()
        self.df_post_flow.plot(color='black',alpha=0.2, ax= ax2)
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Total Flow [mm/hour]')
        ax2.set_title('Total Flow for behavioral runs')
        ax2.get_legend().remove()
        f2.savefig(os.path.join(out_file_path, 'flow_range.png'), dpi=600)





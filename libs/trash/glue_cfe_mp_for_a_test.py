import sys

sys.path.append("../libs/cfe/cfe_py")
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

# Global variables
quantiles = [0.10, 0.5, 0.90]


# Global function
def weighted_quantile(
    values, quantiles, sample_weight=None, values_sorted=False, old_style=False
):
    # Code from https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy/32216049
    """Very close to numpy.percentile, but supports weights.
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
    assert np.all(quantiles >= 0) and np.all(
        quantiles <= 1
    ), "quantiles should be in [0, 1]"

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
class MyGLUE(object):
    def __init__(self, out_path="./", config_path_CFE="", nrun=1, eval_criteria=dict()):
        self.out_path = out_path  # Output folder path
        self.nrun = nrun  # Number of runs
        self.var_names = ["Flow", "Soil Moisture Content"]  # Variables to be analyzed
        self.config_path_CFE = config_path_CFE

        # Evaluation criteria (multi-criteria allowed)
        self.eval_criteria = eval_criteria
        print(f"A number of criterion: {len(eval_criteria)}")
        self.eval_names = []
        for i in range(len(eval_criteria)):
            if eval_criteria[i]["metric"] == "season_transition":
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (d2w_start)'
                )
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (d2w_end)'
                )
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (w2d_start)'
                )
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]} (w2d_end)'
                )
            else:
                self.eval_names.append(
                    f'{eval_criteria[i]["metric"]} on {eval_criteria[i]["variable_to_analyze"]}'
                )
            print(
                f'[{i + 1}] {eval_criteria[i]["metric"]}-based analysis on {eval_criteria[i]["variable_to_analyze"]}'
            )

    def simulation(self, sampled_params):
        # ===============================================================
        # The main model run
        # ===============================================================

        # Preparations
        nth_run = sampled_params[0]
        sampled_params_set = sampled_params[1]
        print(f"Processing {nth_run}/{self.nrun-1}")

        # ===============================================================
        # Write the randomly-generated parameters to the config json file
        # ===============================================================

        # CFE model instance
        template_config_CFE = self.config_path_CFE
        target_config_CFE = os.path.join(
            r"..\2_data_input\Mahurangi\parameters", f"config_cfe_{nth_run}.json"
        )
        shutil.copyfile(template_config_CFE, target_config_CFE)

        # Get the model config file
        with open(target_config_CFE) as data_file:
            self.cfe_cfg = json.load(data_file)

        # Overwrite the model config file
        for i in range(len(sampled_params_set)):
            if sampled_params_set[i][1] in [
                "bb",
                "satdk",
                "slop",
                "satpsi",
                "smcmax",
                "wltsmc",
                "D",
            ]:
                self.cfe_cfg["soil_params"][
                    sampled_params_set[i][1]
                ] = sampled_params_set[i][0]
            else:
                self.cfe_cfg[sampled_params_set[i][1]] = sampled_params_set[i][0]

        # Save the config file with new parameters
        with open(target_config_CFE, "w") as out_file:
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
            sim["Time"] = pd.to_datetime(
                sim["Time"], format="%Y-%m-%d %H:%M:%S"
            )  # Works specifically for CFE

            # Get the comparison data
            obs = obs0[["Time", var_name]].copy()
            obs["Time"] = pd.to_datetime(
                obs["Time"], format="%m/%d/%Y %H:%M"
            )  # Works specifically for Mahurangi data
            # obs["Time"] = pd.to_datetime(obs["Time"], format="%d-%m-%Y %H:%M:%S")

            # Merge observed and simulated timeseries
            df = pd.merge_asof(sim, obs, on="Time")

            sim_synced[var_name] = df[var_name + "_x"].copy()
            obs_synced[var_name] = df[var_name + "_y"].copy()

        # ===============================================================
        # Evalute the outputs
        # ===============================================================

        # Preparation
        eval_result_for_a_run = []
        behavioral_flag = [False] * len(self.eval_criteria)

        # Loop for all evaluation metrics (multi-criteria).
        # Calculate the metrics and judge behavioral vs. non-behavioral
        for i in range(len(self.eval_criteria)):
            # Nash-Sutcliffe scores
            if self.eval_criteria[i]["metric"] == "NSE":
                metric_value = spotpy.objectivefunctions.nashsutcliffe(
                    obs_synced[self.eval_criteria[i]["variable_to_analyze"]],
                    sim_synced[self.eval_criteria[i]["variable_to_analyze"]],
                )
                if metric_value > self.eval_criteria[i]["threshold"]:
                    behavioral_flag[i] = True

            # Kling-Gupta Efficiency scores
            # 　Kling-Gupta efficiencies range from -Inf to 1. Essentially, the closer to 1, the more accurate the model is
            elif self.eval_criteria[i]["metric"] == "KGE":
                metric_value = spotpy.objectivefunctions.kge(
                    obs_synced[self.eval_criteria[i]["variable_to_analyze"]],
                    sim_synced[self.eval_criteria[i]["variable_to_analyze"]],
                )
                if metric_value > self.eval_criteria[i]["threshold"]:
                    behavioral_flag[i] = True

            # Seasonal transition dates
            elif self.eval_criteria[i]["metric"] == "season_transition":
                # Calculate metrics for OBSERVED timeseries as a baseline performance. Run only once
                # if n == 0:
                sig_obs = SMSig(
                    ts_time=df["Time"].to_numpy(),
                    ts_value=obs_synced[
                        self.eval_criteria[i]["variable_to_analyze"]
                    ].to_numpy(),
                    plot_results=False,
                    plot_label="obs",
                )
                # sig_obs.detrend() # TODO:debug
                sig_obs.movmean()
                t_valley = sig_obs.calc_sinecurve()
                season_trans_obs, _, _ = sig_obs.calc_seasontrans(t_valley=t_valley)

                # Calculate metrics for SIMULATED timeseries
                sig_sim = SMSig(
                    ts_time=df["Time"].to_numpy(),
                    ts_value=sim_synced[
                        self.eval_criteria[i]["variable_to_analyze"]
                    ].to_numpy(),
                    plot_results=False,
                    plot_label="sim",
                )
                sig_sim.movmean()
                season_trans_sim, _, _ = sig_sim.calc_seasontrans(t_valley=t_valley)

                # Get the deviations in seasonal transition dates between simulated and observed timeseries
                diff = season_trans_sim - season_trans_obs
                metric_value = abs(np.nanmean(diff, axis=0))
                if all(metric_value < self.eval_criteria[i]["threshold"]):
                    behavioral_flag[i] = True

            # Store evaluation metrics for all criteria for one run
            eval_result_for_a_run.append(metric_value)

        # ===============================================================
        # Judge behavioral vs. non-behavioral
        # ===============================================================

        if all(behavioral_flag):
            # If all criteria is TRUE, the model is behavioral
            result_glue = "Behavioral"
            results = [
                sim_synced["Flow"],
                sim_synced["Soil Moisture Content"],
                nth_run,
                sampled_params_set,
                eval_result_for_a_run,
                result_glue,
            ]
        else:
            # Discard non-behavioral runs
            result_glue = "Non-behavioral"
            results = [
                pd.DataFrame({"Flow": []}),
                pd.DataFrame({"Soil Moisture Content": []}),
                nth_run,
                sampled_params_set,
                eval_result_for_a_run,
                result_glue,
            ]

        print(f"{nth_run}-th run/{self.nrun-1}: {result_glue}")
        print(eval_result_for_a_run)

        return results

    def save_results_to_df(self, all_results):
        # ===============================================================
        # Save results from all runs
        # ===============================================================

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
            sim["Time"] = pd.to_datetime(
                sim["Time"], format="%Y-%m-%d %H:%M:%S"
            )  # Works specifically for CFE

            # Get the comparison data
            obs = obs0[["Time", var_name]].copy()
            obs["Time"] = pd.to_datetime(
                obs["Time"], format="%m/%d/%Y %H:%M"
            )  # Works specifically for Mahurangi data

            # Merge observed and simulated timeseries
            df = pd.merge_asof(sim, obs, on="Time")
            self.df_timeaxis = df["Time"]
            obs_synced[var_name] = df[var_name + "_y"].copy()
        self.obs_synced = obs_synced

        self.df_obs_Q = pd.DataFrame(self.obs_synced["Flow"])
        self.df_obs_Q.set_axis(self.df_timeaxis, axis=0, inplace=True)
        self.df_obs_SM = pd.DataFrame(self.obs_synced["Soil Moisture Content"])
        self.df_obs_SM.set_axis(self.df_timeaxis, axis=0, inplace=True)

        # ==================================================
        # Store results to the dataframe
        # ==================================================

        ## Store GLUE results (behavioral vs. non-behavioral)
        self.glue_results = [np.nan] * len(all_results)
        self.run_id = [np.nan] * len(all_results)
        for i in range(len(all_results)):
            self.run_id[i] = all_results[i][2]
            if all_results[i][5] == "Behavioral":
                boolean_glue_result = True
            else:
                boolean_glue_result = False
            self.glue_results[i] = boolean_glue_result
        self.df_glue_results = pd.DataFrame(
            self.glue_results, index=self.run_id, columns=["Behavioral"]
        )
        behavioral_run_id_index = self.df_glue_results.index[
            self.df_glue_results["Behavioral"].values
        ]

        ## Store timeseries of flow and soil moisture
        # Flow & soil moisture
        for i in range(len(all_results)):
            if self.glue_results[i]:
                if not hasattr(self, "df_behavioral_Q"):
                    self.df_behavioral_Q = all_results[i][0]
                    self.df_behavioral_SM = all_results[i][1]
                else:
                    self.df_behavioral_Q = pd.concat(
                        [self.df_behavioral_Q, all_results[i][0]], axis=1
                    )
                    self.df_behavioral_SM = pd.concat(
                        [self.df_behavioral_SM, all_results[i][1]], axis=1
                    )

        # Store Behavioral runs for Flow and soil moisture
        self.run_id_behavioral = [
            self.run_id[i] for i in range(len(self.run_id)) if self.glue_results[i]
        ]
        self.df_behavioral_Q.set_axis(self.run_id_behavioral, axis=1, inplace=True)
        self.df_behavioral_Q.set_axis(self.df_timeaxis, axis=0, inplace=True)
        self.df_behavioral_SM.set_axis(self.run_id_behavioral, axis=1, inplace=True)
        self.df_behavioral_SM.set_axis(self.df_timeaxis, axis=0, inplace=True)

        ## Store PRIOR parameters
        self.pri_paras = [np.nan] * len(all_results)
        for i in range(len(all_results)):
            self.pri_paras[i] = all_results[i][3]

        param_names = []
        for i in range(len(self.pri_paras[0])):
            param_names.append(self.pri_paras[0][i][1])
        param_values = np.empty((self.nrun, len(param_names)))
        param_values[:] = np.nan
        for i in range(len(param_names)):
            for j in range(self.nrun):
                param_values[j][i] = self.pri_paras[j][i][0]
        self.df_pri_paras = pd.DataFrame(param_values, columns=param_names)

        # Store POSTERIOR paramters for behavioral runs
        self.df_post_paras = self.df_pri_paras.iloc[behavioral_run_id_index].copy()

        ## Store Evaluation metrics for all runs
        self.eval = [np.nan] * len(all_results)
        for i in range(len(all_results)):
            self.eval[i] = all_results[i][4]

        eval_values = np.empty((len(all_results), len(self.eval_names)))
        eval_values[:] = np.nan
        for j in range(len(self.eval)):
            season_idx = -9999
            for i in range(len(self.eval_names)):
                if "season_transition" in self.eval_names[i]:
                    if season_idx == -9999:
                        season_idx = i
                        season_num = 0
                    eval_values[j][i] = self.eval[j][season_idx][season_num]
                    season_num += 1
                else:
                    eval_values[j][i] = self.eval[j][i]

        # Store behavioral evaluations
        self.df_eval = pd.DataFrame(
            eval_values, index=self.run_id, columns=self.eval_names
        )
        self.df_post_eval = self.df_eval.iloc[behavioral_run_id_index].copy()

    def post_process(self):
        # post-process the results
        print("--- Post-processing the simulated results ---")

        # Calculate weights
        weight = np.empty((len(self.df_post_eval), len(self.eval_names)))
        j = int(0)
        # Loop for all evaluation metrics
        for i in range(len(self.eval_criteria)):
            if (
                self.eval_criteria[i]["metric"] == "NSE"
                or self.eval_criteria[i]["metric"] == "KGE"
            ):
                # For Nash-Sutcliffe and Kling-Gupta Efficiency scores
                weight[:, j] = (
                    (
                        self.df_post_eval[self.eval_names[j]]
                        - self.eval_criteria[i]["threshold"]
                    )
                    / sum(
                        self.df_post_eval[self.eval_names[j]]
                        - self.eval_criteria[i]["threshold"]
                    )
                ).to_numpy()
                j += int(1)
            elif self.eval_criteria[i]["metric"] == "season_transition":
                for k in range(4):
                    # For seasonal transition dates
                    weight[:, j + k] = triangle_weight(
                        self.df_post_eval[self.eval_names[j + k]],
                        a=-1 * self.eval_criteria[i]["threshold"],
                        b=0,
                        c=self.eval_criteria[i]["threshold"],
                    )
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
                quantile[t, :] = weighted_quantile(
                    values=values,
                    quantiles=quantiles,
                    sample_weight=avg_weight,
                    values_sorted=False,
                    old_style=False,
                )
            df_simrange = pd.DataFrame(
                quantile,
                index=df_behavioral.index,
                columns=["lowerlim", "median", "upperlim"],
            )
            if var_name == "Flow":
                self.df_Q_simrange = df_simrange.copy()
            elif var_name == "Soil Moisture Content":
                self.df_SM_simrange = df_simrange.copy()

    def to_csv(self, df_param_to_calibrate=None):
        print("--- Saving data into csv file ---")

        df_param_to_calibrate.to_csv(
            os.path.join(self.out_path, "parameter_bounds_used.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )
        self.df_glue_results.to_csv(
            os.path.join(self.out_path, "glue_results.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )
        self.df_pri_paras.to_csv(
            os.path.join(self.out_path, "paramter_priori.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )
        self.df_eval.to_csv(
            os.path.join(self.out_path, "evaluations.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )
        if hasattr(self, "df_behavioral_Q"):
            self.df_behavioral_Q.to_csv(
                os.path.join(self.out_path, "behavioral_Q.csv"),
                sep=",",
                header=True,
                index=True,
                encoding="utf-8",
                na_rep="nan",
            )
        if hasattr(self, "df_behavioral_SM"):
            self.df_behavioral_SM.to_csv(
                os.path.join(self.out_path, "behavioral_SM.csv"),
                sep=",",
                header=True,
                index=True,
                encoding="utf-8",
                na_rep="nan",
            )
        if hasattr(self, "df_Q_simrange"):
            self.df_Q_simrange.to_csv(
                os.path.join(self.out_path, "quantiles_Q.csv"),
                sep=",",
                header=True,
                index=True,
                encoding="utf-8",
                na_rep="nan",
            )
        if hasattr(self, "df_SM_simrange"):
            self.df_SM_simrange.to_csv(
                os.path.join(self.out_path, "quantiles_SM.csv"),
                sep=",",
                header=True,
                index=True,
                encoding="utf-8",
                na_rep="nan",
            )

    def plot(self, plot_type=None):
        # Plot the results
        print("--- Saving the plots ---")

        # Histogram
        # Prior vs. posterior parameter distributions
        if plot_type == "param_hist":
            nparas = len(self.df_pri_paras.columns)
            f = plt.figure(figsize=(4 * 4, 4 * 4))

            for i in range(nparas):
                target_para = self.df_pri_paras.columns[i]
                ax1 = f.add_subplot(4, 4, i + 1)

                self.df_pri_paras[target_para].plot.hist(
                    bins=10, alpha=0.4, ax=ax1, color="#3182bd", label="Prior"
                )

                if hasattr(self, "df_post_paras"):
                    self.df_post_paras[target_para].plot.hist(
                        bins=10, alpha=0.8, ax=ax1, color="#3182bd", label="Posterior"
                    )

                ax1.set_xlabel(target_para)
                ax1.legend()

                if i != 0:
                    ax1.yaxis.set_visible(False)

            # f.plot()

            f.savefig(os.path.join(self.out_path, "param_dist.png"), dpi=600)

        # Dotty plot
        # Prior vs. posterior parameter distributions
        if plot_type == "dotty":
            if hasattr(self, "df_post_paras"):
                nparas = len(self.df_pri_paras.columns)

                for j in range(len(self.df_post_eval.columns)):
                    f = plt.figure(figsize=(4 * 4, 4 * 4), constrained_layout=True)
                    f.tight_layout()
                    target_eval = self.df_post_eval.columns[j]
                    for i in range(nparas):
                        target_para = self.df_pri_paras.columns[i]
                        ax1 = f.add_subplot(4, 4, i + 1)
                        ax1.scatter(
                            self.df_post_paras[target_para],
                            self.df_post_eval[target_eval],
                            alpha=0.5,
                        )
                        ax1.tick_params(axis="both", which="major", labelsize=16)
                        ax1.set_xlabel(target_para, fontsize=16)
                        ax1.set_ylabel(target_eval, fontsize=16)

                    # if i !=0:
                    #     ax1.yaxis.set_visible(False)

                    # f.plot()
                    f.savefig(
                        os.path.join(
                            self.out_path, "param_dotty_%s.png" % (target_eval)
                        ),
                        dpi=600,
                    )

        # Dotty plot
        # Parameter interactions for the behavioral runs
        if plot_type == "dotty_interaction":
            param_interset = [
                "bb",
                "satdk",
                "slop",
                "satpsi",
                "smcmax",
                "wltsmc",
                "alpha_fc",
                "lksatfac",
                "D",
                "trigger_z_fact",
                "max_gw_storage",
                "Cgw",
                "expon",
                "refkdt",
                "K_nash",
            ]

            if hasattr(self, "df_post_paras"):
                f = plt.figure(figsize=(20, 20), constrained_layout=True)
                f.tight_layout()
                n_plot = 0
                for i in range(len(param_interset)):
                    for j in range(len(param_interset)):
                        n_plot += 1
                        para0 = param_interset[j]
                        para1 = param_interset[i]
                        ax1 = f.add_subplot(
                            len(param_interset), len(param_interset), n_plot
                        )
                        x = self.df_post_paras[para0]
                        y = self.df_post_paras[para1]
                        ax1.scatter(x, y, alpha=0.5)

                        if i == 0:
                            ax1.xaxis.set_label_position("top")
                            ax1.set_xlabel(para0)
                        if j == 0:
                            ax1.set_ylabel(para1)
                        ax1.tick_params(direction="in")

                f.savefig(
                    os.path.join(self.out_path, "param_dotty_interaction.png"), dpi=600
                )

        # Time series of data
        # Flow and Soil Moisture Content
        if plot_type == "timeseries":
            for var_name in self.var_names:
                if var_name == "Flow":
                    df_simrange = self.df_Q_simrange
                    df_obs = self.df_obs_Q
                    obs_label = "Observed flow"
                    y_label = "Total flow [mm/hour]"
                    title = "Total flow for behavioral runs"
                    fn = "Q_range.png"
                elif var_name == "Soil Moisture Content":
                    var_name = "Soil Moisture Content"
                    df_simrange = self.df_SM_simrange
                    df_obs = self.df_obs_SM
                    obs_label = "Observed soil moisture"
                    y_label = "Volumetric Soil Moisture Content [m^3/m^3]"
                    title = "Soil moisture for behavioral runs"
                    fn = "SM_range.png"

                f2 = plt.figure(figsize=(8, 6))
                ax2 = f2.add_subplot()

                df_simrange["lowerlim"].plot(
                    color="black",
                    alpha=0.2,
                    ax=ax2,
                    label=f"{quantiles[0]*100} percentile",
                )
                df_simrange["upperlim"].plot(
                    color="black",
                    alpha=0.2,
                    ax=ax2,
                    label=f"{quantiles[2]*100} percentile",
                )
                df_obs[var_name].plot(color="black", alpha=1, ax=ax2, label=obs_label)
                plt.fill_between(
                    df_simrange.index,
                    df_simrange["upperlim"],
                    df_simrange["lowerlim"],
                    facecolor="green",
                    alpha=0.2,
                    interpolate=True,
                    label="Predicted range",
                )
                if var_name == "Flow":
                    ax2.set_yscale("log")
                    ax2.set_ylim([1e-07, 1e-01])
                ax2.set_xlabel("Time")
                ax2.set_ylabel(y_label)
                ax2.set_title(title)
                ax2.legend()
                f2.autofmt_xdate()
                f2.savefig(os.path.join(self.out_path, fn), dpi=600)

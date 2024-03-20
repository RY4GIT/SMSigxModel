import os
import sys

sys.path.append(os.path.join(os.getcwd(), "libs", "SMSig"))
from sig_seasontrans import SMSig

import numpy as np
import pandas as pd
from statistics import median
import json
import spotpy


class Evaluator:
    def __init__(self, config=None, simulation=None, observation=None) -> None:
        self.flow_var_name = "Flow"
        self.soilmoisture_var_name = "Soil Moisture Content"

        self.config = config
        self.observation = observation
        self.simulation = simulation

    def evaluate(self, nth_run=None):
        """Evaluate simulation results."""

        eval_hourly = self._evaluate_hourly()
        eval_monthly = self._evaluate_monthly()

        # Convert to desired format
        eval_hourly["run_id"] = nth_run
        eval_hourly.set_index("run_id", inplace=True)

        eval_monthly.reset_index(inplace=True)
        eval_monthly["run_id"] = nth_run
        eval_monthly.set_index("run_id", inplace=True)

        return eval_hourly, eval_monthly

    def _evaluate_hourly(self):
        """Evaluate metrics at model (hourly) timestep."""

        metrics = {
            "NSE on Flow": self.calc_NSE(self.flow_var_name),
            "NSE on Soil": self.calc_NSE(self.soilmoisture_var_name),
            "KGE on Flow": self.calc_KGE(self.flow_var_name),
            "KGE on Soil": self.calc_KGE(self.soilmoisture_var_name),
        }

        diff_season_trans = self.calc_SeasonTrans(self.soilmoisture_var_name)
        abs_diff_SeasonTransDate = abs(np.nanmean(diff_season_trans, axis=0))
        metrics.update(
            {
                "SeasonTrans of Soil dry2wet_start": abs_diff_SeasonTransDate[0],
                "SeasonTrans of Soil dry2wet_end": abs_diff_SeasonTransDate[1],
                "SeasonTrans of Soil wet2dry_start": abs_diff_SeasonTransDate[2],
                "SeasonTrans of Soil wet2dry_end": abs_diff_SeasonTransDate[3],
                "SeasonTrans of Soil dry2wet_start_raw": [
                    diff_season_trans[:, 0].tolist()
                ],
                "SeasonTrans of Soil dry2wet_end_raw": [
                    diff_season_trans[:, 1].tolist()
                ],
                "SeasonTrans of Soil wet2dry_start_raw": [
                    diff_season_trans[:, 2].tolist()
                ],
                "SeasonTrans of Soil wet2dry_end_raw": [
                    diff_season_trans[:, 3].tolist()
                ],
            }
        )

        return pd.DataFrame(metrics, index=[0])

    def _evaluate_monthly(self):
        """Evaluate metrics at model (monthly) timestep."""

        self.observation_monthly = self.df_to_monthly_timestep(self.observation)
        self.simulation_monthly = self.df_to_monthly_timestep(self.simulation)

        median_flow_obs = median(self.observation["Flow"])
        high_flow_obs = 9 * median_flow_obs

        metrics = {
            "Q_mean_obs": self.calc_monthly_Q_mean(
                self.flow_var_name, self.observation_monthly
            ),
            "Q_mean_sim": self.calc_monthly_Q_mean(
                self.flow_var_name, self.simulation_monthly
            ),
            "high_flow_freq_obs": self.calc_monthly_high_flow_freq(
                self.flow_var_name, self.observation, high_flow_obs
            ),
            "high_flow_freq_sim": self.calc_monthly_high_flow_freq(
                self.flow_var_name, self.simulation, high_flow_obs
            ),
            "RR_obs": self.calc_monthly_RR(
                self.flow_var_name,
                self.observation_monthly,
                self.observation["Rainfall"],
            ),
            "RR_sim": self.calc_monthly_RR(
                self.flow_var_name,
                self.simulation_monthly,
                self.observation["Rainfall"],
            ),
        }

        return pd.concat(metrics.values(), axis=1, keys=metrics.keys())

    def calc_NSE(self, variable):
        obs = self.observation[variable]
        sim = self.simulation[variable]
        return spotpy.objectivefunctions.nashsutcliffe(
            obs[~np.isnan(obs)], sim[~np.isnan(obs)]
        )

    def calc_KGE(self, variable):
        obs = self.observation[variable]
        sim = self.simulation[variable]
        return spotpy.objectivefunctions.kge(obs[~np.isnan(obs)], sim[~np.isnan(obs)])

    def calc_SeasonTrans(self, variable):
        # Calculate metrics for observed timeseries
        sig_obs = SMSig(
            t=self.observation.index.to_numpy(),
            sm=self.observation[variable].to_numpy(),
            plot_results=False,
            verbose=False,
        )

        seasonal_cycle = self.get_seasonal_cycle()
        parameter_config = self.get_SeasonTrans_config()
        season_trans_obs = sig_obs.calc_seasontrans(
            seasonal_cycle=seasonal_cycle, parameter_config=parameter_config
        )

        # Calculate metrics for SIMULATED timeseries
        sig_sim = SMSig(
            t=self.simulation.index.to_numpy(),
            sm=self.simulation[variable].to_numpy(),
            plot_results=False,
            verbose=False,
        )
        season_trans_sim = sig_sim.calc_seasontrans(
            seasonal_cycle=seasonal_cycle, parameter_config=parameter_config
        )

        # Get the deviations in seasonal transition dates between simulated and observed timeseries
        diff = season_trans_sim - season_trans_obs

        return diff

    def get_seasonal_cycle(self):
        data_dir = os.path.dirname(self.config["PATHS"]["cfe_config"])
        seasonal_cycle = pd.read_csv(
            os.path.join(data_dir, "seasonal_cycle.csv"),
            parse_dates=["start_date", "end_date"],
        )
        return seasonal_cycle

    def get_SeasonTrans_config(self):
        config_path = os.path.join(
            os.path.dirname(self.config["PATHS"]["cfe_config"]),
            "seasonal_transition_config.json",
        )
        with open(config_path, "r") as config_file:
            config = json.load(config_file)
        return config

    def df_to_monthly_timestep(self, df):
        df_monthly = df.resample("M").sum().copy()
        df_monthly.drop(df_monthly.tail(1).index, inplace=True)
        df_monthly.drop(df_monthly.head(1).index, inplace=True)
        return df_monthly

    def calc_monthly_Q_mean(self, variable, df):
        return df[variable].resample("M").mean()

    def calc_monthly_high_flow_freq(self, variable, df, high_flow_criteria):
        monthly_high_flow_freq = (
            df[variable]
            .resample("M")
            .apply(lambda x: sum(x > high_flow_criteria) / len(x))
        )
        return monthly_high_flow_freq

    def calc_monthly_RR(self, variable, df, precip):
        monthly_flow = df[variable].resample("M").sum()
        monthly_precip = precip.resample("M").sum()
        return monthly_flow / monthly_precip

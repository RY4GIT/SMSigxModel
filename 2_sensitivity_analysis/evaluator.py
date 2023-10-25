import os
import sys

sys.path.append(os.path.join(os.getcwd(), "libs", "SMSig"))
from sig_seasontrans import SMSig

import numpy as np
import pandas as pd

import spotpy


class Evaluator:
    def __init__(
        self, config=None, simulation=None, observation=None, time_index=None
    ) -> None:
        self.config = config
        self.observation = observation
        self.simulation = simulation
        self.time_index = time_index

        self.like_measure = self.config["SALib"]["like_measure"]

    def evaluate(self):
        """Evaluate simulation results."""

        # Calculate objective metrics
        if self.like_measure == "NashSutcliffe":
            like = self.calc_NSE()
        elif self.like_measure == "KGE":
            like = self.calc_KGE()
        elif "SeasonTrans" in self.like_measure:
            season_trans = self.calc_SeasonTrans()

            if self.like_measure == "SeasonTrans of Soil dry2wet_start":
                like = season_trans[0]
            elif self.like_measure == "SeasonTrans of Soil dry2wet_end":
                like = season_trans[1]
            elif self.like_measure == "SeasonTrans of Soil wet2dry_start":
                like = season_trans[2]
            elif self.like_measure == "SeasonTrans of Soil wet2dry_end":
                like = season_trans[3]

        return like

    def calc_NSE(self):
        obs = self.observation
        sim = self.simulation
        return spotpy.objectivefunctions.nashsutcliffe(
            obs[~np.isnan(obs)], sim[~np.isnan(obs)]
        )

    def calc_KGE(self):
        obs = self.observation
        sim = self.simulation
        return spotpy.objectivefunctions.kge(obs[~np.isnan(obs)], sim[~np.isnan(obs)])

    def calc_SeasonTrans(self):
        # Calculate metrics for observed timeseries
        sig_obs = SMSig(
            ts_time=self.time_index,
            ts_value=self.observation.to_numpy(),
            plot_results=False,
            plot_label="obs",
        )

        t_valley = self.get_t_valley()
        season_trans_obs, _, _ = sig_obs.calc_seasontrans(t_valley=t_valley)

        # Calculate metrics for SIMULATED timeseries
        sig_sim = SMSig(
            ts_time=self.time_index,
            ts_value=self.simulation.to_numpy(),
            plot_results=False,
            plot_label="sim",
        )
        season_trans_sim, _, _ = sig_sim.calc_seasontrans(t_valley=t_valley)

        # Get the deviations in seasonal transition dates between simulated and observed timeseries
        diff = season_trans_sim - season_trans_obs
        abs_diff_SeasonTransDate = abs(np.nanmean(diff, axis=0))

        return abs_diff_SeasonTransDate

    def get_t_valley(self):
        data_directory = os.path.dirname(self.config["PATHS"]["cfe_config"])
        _t_valley_manual_input = pd.read_csv(
            os.path.join(data_directory, "seasonal_cycel_valleys.csv"), header=None
        )
        t_valley_manual_input = pd.to_datetime(_t_valley_manual_input[0])
        return t_valley_manual_input

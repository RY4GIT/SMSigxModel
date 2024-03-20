import os
import sys

sys.path.append(os.path.join(os.getcwd(), "libs", "SMSig"))
from sig_seasontrans import SMSig

import numpy as np
import pandas as pd
import json
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
            t=self.time_index,
            sm=self.observation.to_numpy(),
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
            t=self.time_index,
            sm=self.simulation.to_numpy(),
            plot_results=False,
            verbose=False,
        )
        season_trans_sim = sig_sim.calc_seasontrans(
            seasonal_cycle=seasonal_cycle, parameter_config=parameter_config
        )

        # Get the deviations in seasonal transition dates between simulated and observed timeseries
        diff = season_trans_sim - season_trans_obs
        abs_diff_SeasonTransDate = abs(np.nanmean(diff, axis=0))

        return abs_diff_SeasonTransDate

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

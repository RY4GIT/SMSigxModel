import os

from model import CFEmodel
from evaluator import Evaluator

import spotpy
import os
import pandas as pd
import numpy as np
import datetime
from sampler import mylhs, spot_setup


# GLUE object
class Agent_GLUE_CFE(object):
    def __init__(self, config=None):
        self.config = config

        # Define output folder
        # Get the current date in YYYY-MM-DD format
        current_date = datetime.date.today().strftime("%Y-%m-%d")
        self.out_path = os.path.join(
            config["PATHS"]["homedir"],
            "results",
            f"{self.config['DATA']['site']}-{current_date}",
        )
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)

        # Define the random seed for reproducibility
        np.random.seed(0)

        # Define the GLUE number of runs
        if config["Multiprocessing"]["run_multithreadding"] == "False":
            self.nrun = 1
        elif config["Multiprocessing"]["run_multithreadding"] == "True":
            self.nrun = int(config["GLUE"]["nrun"])

    def generate_params(self):
        """Generate random parameters for GLUE run

        Returns:
            sampled_params: randomly sampled parameters
        """

        # Read parameter bounds defined from an excel file
        df = pd.read_csv(self.config["PATHS"]["GLUE_config"])
        self.df_param_to_calibrate = df[df["calibrate"] == 1]

        # Generate random parameters from uniform distribution

        print("--- Sampling parameters ---")
        # Latin Hypercube Sampling (LHS) apprach
        spot_setup_instance = spot_setup(
            df_param_to_calibrate=self.df_param_to_calibrate
        )
        sampler = mylhs(spot_setup_instance)
        sampled_params = sampler.sample(self.nrun)

        # ######### OPTION ##########
        # # Monte Carlo Sampling (MCRS) approach
        # sampled_params = [
        #     [i, spotpy.parameter.generate(spot_setup_instance.params)]
        #     for i in range(self.nrun)
        # ]
        # ###########################

        # alternatively, use Latin Hypercube Sampling (LHS)
        # https://github.com/thouska/spotpy/blob/master/src/spotpy/algorithms/lhs.py
        # Store parameter bounds used in the analysis

        self.df_param_to_calibrate.to_csv(
            os.path.join(self.out_path, "parameter_bounds_used.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )

        return sampled_params

    def simulation(self, sampled_param_set):
        """One CFE run and evaluation for a sampled parameter set"""

        # Preparations
        nth_run = sampled_param_set[0]
        sampled_params_set = sampled_param_set[1:3]
        print(f"Processing {nth_run}/{self.nrun-1}")

        # Run CFE
        cfe_instance = CFEmodel(config=self.config)
        cfe_instance.initialize(nth_run=nth_run, sampled_params_set=sampled_params_set)
        obs, sim = cfe_instance.run()

        # Evaluate the outputs
        evaluator = Evaluator(config=self.config, observation=obs, simulation=sim)
        eval_metrics, eval_metrics_monthly = evaluator.evaluate(nth_run=nth_run)

        print(f"Evaluation {nth_run}/{self.nrun-1}")
        print(eval_metrics.to_string(index=False))

        results = [nth_run, sampled_params_set, eval_metrics, eval_metrics_monthly]

        return results

    def finalize(self, all_results=None):
        """Join and save results from all runs to dataframe"""

        # Store run ID
        self.run_id = [result[0] for result in all_results]

        # Store prior parameters
        self._pri_paras = [result[1] for result in all_results]
        param_names = [param[1] for param in self._pri_paras[0]]
        param_values = np.array(
            [
                [param_set[i][0] for i in range(len(param_names))]
                for param_set in self._pri_paras
            ]
        )
        self.prior_params = pd.DataFrame(param_values, columns=param_names)
        self.prior_params["run_id"] = self.run_id
        self.prior_params.set_index("run_id", inplace=True)
        self.prior_params.to_csv(
            os.path.join(self.out_path, "prior_parameters.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )

        # Store Evaluation metrics (whole period)
        self.eval_metrics = pd.concat([result[2] for result in all_results])
        self.eval_metrics.to_csv(
            os.path.join(self.out_path, "evaluation_metrics.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )

        # Store Evaluation metrics (monthly)
        self.eval_metrics_mo = pd.concat([result[3] for result in all_results])
        self.eval_metrics_mo.to_csv(
            os.path.join(self.out_path, "evaluation_metrics_monthly.csv"),
            sep=",",
            header=True,
            index=True,
            encoding="utf-8",
            na_rep="nan",
        )

        print(f"--- Saved results to {self.out_path} ---")

import os
import spotpy
import pandas as pd

pd.options.mode.chained_assignment = None
import json
import warnings

# %matplotlib inline
import itertools
from math import pi
from matplotlib.legend_handler import HandlerPatch

import sys

sys.path.append(os.path.join(os.getcwd(), "libs", "cfe_py"))
from bmi_cfe import BMI_CFE

import shutil
import random

from evaluator import Evaluator


def duplicate_file(source_path):
    i = random.randint(1, 9999)
    # Determine the directory and make a path for the new file
    directory = os.path.join(
        os.path.dirname(source_path),
        "temporary_parameter_files_for_sensitivity_analysis",
    )
    destination_path = os.path.join(directory, f"cfe_config_{i}.json")

    # Copy the source file to the new location
    shutil.copy2(source_path, destination_path)

    return destination_path  # Optionally return the path to the new file


class CFEmodel:
    def __init__(self, config=None, problem=None, X=None):
        # Configs
        self.config = config
        self.like_measure = config["SALib"]["like_measure"]
        self.eval_variable = config["SALib"]["eval_variable"]
        self.warmup_offset = int(config["CFE"]["warmup_offset"])
        self.warmup_iteration = int(config["CFE"]["warmup_iteration"])

        # Copy the CFE config file for sensitivty analysis
        destination_path = duplicate_file(
            source_path=self.config["PATHS"]["cfe_config"]
        )

        self.cfe_instance = BMI_CFE(cfg_file=destination_path)

        # write the randomly-generated parameters to the config json file
        with open(self.cfe_instance.cfg_file) as data_file:
            cfe_cfg = json.load(data_file)

        for i, param_name in enumerate(problem["names"]):
            if param_name in [
                "bb",
                "satdk",
                "slop",
                "satpsi",
                "smcmax",
                "wltsmc",
                "D",
            ]:
                cfe_cfg["soil_params"][param_name] = X[i]
            else:
                cfe_cfg[param_name] = X[i]

        with open(self.cfe_instance.cfg_file, "w") as out_file:
            json.dump(cfe_cfg, out_file, indent=4)

        # Here the model is actualy started with a unique parameter combination that it gets from spotpy for each time the model is called
        self.cfe_instance.initialize()

    def run(self):
        self.cfe_instance.run_unit_test(
            plot=False,
            print_fluxes=False,
            warm_up=True,
            warmup_offset=self.warmup_offset,
            warmup_iteration=self.warmup_iteration,
        )

    def evaluate(self):
        sim = self.to_datetime(
            df=self.cfe_instance.cfe_output_data[["Time", self.eval_variable]],
            time_column="Time",
            format="%Y-%m-%d",
        )
        obs = self.to_datetime(
            df=self.cfe_instance.unit_test_data[["Time", self.eval_variable]],
            time_column="Time",
            format=self.config["DATA"]["evaldata_time_format"],
        )

        if obs.index[0] != sim.index[0]:
            warnings.warn(
                "The start of observation and simulation time is different by %s"
                % obs.index[0]
                - sim.index[0]
            )

        if obs.index[-1] != sim.index[-1]:
            warnings.warn(
                "The end of observation and simulation time is different by %s"
                % obs.index[-1]
                - sim.index[-1]
            )

        df = pd.merge_asof(sim, obs, on="Time")

        sim_synced = df[self.eval_variable + "_x"]
        obs_synced = df[self.eval_variable + "_y"]

        evaluator = Evaluator(
            config=self.config,
            observation=obs_synced,
            simulation=sim_synced,
            time_index=pd.to_datetime(df["Time"]).values,
        )

        like = evaluator.evaluate()
        return like

    def to_datetime(self, df, time_column, format):
        df = df.copy()
        df[time_column] = pd.to_datetime(df[time_column], format=format)
        return df.set_index(time_column)

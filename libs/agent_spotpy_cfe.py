import spotpy
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
from .py_cfe import cfe
import os
import pandas as pd
import numpy as np
import json

class spot_setup(object):
    def __init__(self, cfe_input, obj_func=None):

        # Find Path to CFE on users system
        self.obj_func = False

        # Model
        self.myCFE = cfe_input

        # define parameter bounds
        self.params = [
            spotpy.parameter.Uniform('maxsmc', low=0.1, high=0.7),
            spotpy.parameter.Uniform('wltsmc', low=0.0, high=0.6)
            ]
        #     spotpy.parameter.Uniform('satdk', low=1e-07, high=1e-05)
        # ]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, x):
        # write the randomly-generated parameters to the config json file
        with open(self.myCFE.cfg_file) as data_file:
            self.cfe_cfg = json.load(data_file)

        self.cfe_cfg["soil_params"]['maxsmc'] = x['maxsmc']
        self.cfe_cfg["soil_params"]['wltsmc'] = x['wltsmc']
        # self.cfe_cfg["soil_params"]['satdk'] = x['satdk']

        with open(self.myCFE.cfg_file, 'w') as out_file:
            json.dump(self.cfe_cfg, out_file)

        # Here the model is actualy started with a unique parameter combination that it gets from spotpy for each time the model is called
        data = self.myCFE.run_unit_test(plot=False)
        sim = data[["Time", "Total Discharge"]]
        sim["Time"] = pd.to_datetime(sim["Time"])
        return sim

    def evaluation(self):
        # Get the comparison data
        data = self.myCFE.load_unit_test_data()
        obs = data[["Time", "Total Discharge"]]
        obs["Time"] = pd.to_datetime(obs["Time"])
        return obs

    def objectivefunction(self, simulation, evaluation, params=None):
        df = pd.merge_asof(simulation, evaluation, on = "Time")
        sim_synced = df["Total Discharge_x"]
        obs_synced = df["Total Discharge_y"]

        # SPOTPY expects to get one or multiple values back,
        # that define the performance of the model run
        if not self.obj_func:
            # This is used if not overwritten by user
            like = spotpy.objectivefunctions.kge_non_parametric(obs_synced, sim_synced)
        else:
            # Way to ensure flexible spot setup class
            like = self.obj_func(obs_synced, sim_synced)
        return like




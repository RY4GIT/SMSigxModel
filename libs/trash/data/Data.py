
"""A file to store the function where we read the input data"""

import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime


class Data(config=None):
    def __init__(self): 
        
        self.start_time = None
        self.end_time = None
        self.observation = None
        self.forcing = None
        
    def get_forcings(self):
        
    def get_observations(self):
        
    def sync_timeseries(self, observation, simulation):
    
    


        sim = self.cfe_instance.cfe_output_data[["Time", var_measure]]
        sim.loc[:, "Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S")
        sim = sim.set_index("Time")

        # Get the comparison data
        data = self.cfe_instance.unit_test_data
        obs = data[["Time", var_measure]]
        obs.loc[:, "Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M") # "%d-%b-%Y %H:%M:%S"
        obs = obs.set_index("Time")

        if obs.index[0] != sim.index[0]:
            diff_time = obs.index[0] - sim.index[0]
            warnings.warn("The start of observation and simulation time is different by %s" % diff_time)

        if obs.index[-1] != sim.index[-1]:
            diff_time = obs.index[-1] - sim.index[-1]
            warnings.warn("The end of observation and simulation time is different by %s" % diff_time)

        df = pd.merge_asof(sim, obs, on = "Time")

        sim_synced = df[var_measure+"_x"]
        obs_synced = df[var_measure+"_y"]

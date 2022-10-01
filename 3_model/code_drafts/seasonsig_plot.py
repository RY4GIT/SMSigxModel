import spotpy
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import os
import sys
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/cfe/py_cfe")
import cfe
import bmi_cfe

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/SMSig")
from sig_seasontrans import SMSig

# specify current directory create output directory if it does not exist
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model")
os.getcwd()
data_file_path = '..\\2_data_input\\Mahurangi\\full'

myCFE = bmi_cfe.BMI_CFE(os.path.join(data_file_path, 'config_cfe.json'))


"""
# define parameter bounds
params = [
    spotpy.parameter.Uniform('bb', low=2, high=15),
    spotpy.parameter.Uniform('satdk', low=0, high=1),
    spotpy.parameter.Uniform('slop', low=0, high=1),
    spotpy.parameter.Uniform('satpsi', low=0.02, high=0.78),
    spotpy.parameter.Uniform('smcmax', low=0.33, high=0.7),
    spotpy.parameter.Uniform('wltsmc', low=0.0, high=0.57),
    spotpy.parameter.Uniform('exponent_secondary', low=1, high=8),
    spotpy.parameter.Uniform('coeff_secondary', low=0.01, high=3),
    spotpy.parameter.Uniform('trigger_z_m', low=0.01, high=0.87),
    spotpy.parameter.Uniform('fc_atm_press_fraction', low=0.10, high=0.33),
    spotpy.parameter.Uniform('max_gw_storage', low=10, high=250),
    spotpy.parameter.Uniform('Cgw', low=0.01, high=1),
    spotpy.parameter.Uniform('expon', low=1, high=8),
    spotpy.parameter.Uniform('refkdt', low=0.1, high=4)
]


sampled = spotpy.parameter.generate(self.params)

# Get the model config file
with open(self.myCFE.cfg_file) as data_file:
    self.cfe_cfg = json.load(data_file)

# Overwrite the model config file
for i in range(len(self.sampled)):
    if self.sampled[i][1] in ['bb', 'satdk', 'satpsi', 'slop', 'smcmax', 'wltsmc', 'exponent_secondary', 'coeff_secondary', 'D']:
        self.cfe_cfg["soil_params"][self.sampled[i][1]] = self.sampled[i][0]
    else:
        self.cfe_cfg[self.sampled[i][1]] = self.sampled[i][0]

# Save the config file with new parameters
with open(self.myCFE.cfg_file, 'w') as out_file:
    json.dump(self.cfe_cfg, out_file)
"""
# ===============================================================
# Actual model run
# ===============================================================
myCFE.initialize()
sim0 = myCFE.run_unit_test(plot=False)
obs0 = myCFE.load_unit_test_data()
var_name = "Soil Moisture Content"

# Get the simulated data
sim = sim0[["Time", var_name]]
sim["Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S") # Works specifically for CFE

# Get the comparison data
obs = obs0[["Time", var_name]]
obs["Time"] = pd.to_datetime(obs["Time"], format="%d-%b-%Y %H:%M:%S") # Works specifically for Mahurangi data

# Merge observed and simulated timeseries
df = pd.merge_asof(sim, obs, on="Time")
sim_synced = df[var_name + "_x"]
obs_synced = df[var_name + "_y"]


sig_obs = SMSig(ts_time=df["Time"].to_numpy(), ts_value=obs_synced.to_numpy(), plot_results=True, plot_label="obs")
sig_obs.movmean()
t_valley = sig_obs.calc_sinecurve()
season_trans_obs = sig_obs.calc_seasontrans(t_valley=t_valley)

# simulated
sig_sim = SMSig(ts_time=df["Time"].to_numpy(), ts_value=sim_synced.to_numpy(), plot_results=True, plot_label="sim")
sig_sim.movmean()
season_trans_sim = sig_sim.calc_seasontrans(t_valley=t_valley)

diff = season_trans_sim - season_trans_obs
print(diff)





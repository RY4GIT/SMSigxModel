# %% [markdown]
# ## Import libraries

# %%
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np
from datetime import datetime, timedelta

sys.path.append(f"./libs/cfe_py/")
from bmi_cfe import BMI_CFE

from sig_seasontrans import SMSig


# %%
def to_datetime(df, time_column, format="%Y-%m-%d %H:%M:%S"):
    df = df.copy()
    df[time_column] = pd.to_datetime(df[time_column], format=format)
    return df.set_index(time_column)


# %%
home_dir = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis"
data_dir = "data"
site = "Coweeta"
test_file = "test_daily_2014_2018_sm_basinavg.csv"

# %%
_data = pd.read_csv(os.path.join(home_dir, "data", site, test_file))
_data = to_datetime(_data, "Time")
data = _data["Soil Moisture Content"]

data.head()

sig_obs = SMSig(
    t=data.index.to_numpy(), sm=data.to_numpy(), plot_results=True, verbose=True
)

seasonal_cycle = pd.read_csv(
    os.path.join(home_dir, "data", site, "seasonal_cycle.csv"),
    parse_dates=["start_date", "end_date"],
)
_parameter_config = os.path.join(
    home_dir, "data", site, "seasonal_transition_config.json"
)
with open(_parameter_config, "r") as config_file:
    config = json.load(config_file)

season_trans_obs = sig_obs.calc_seasontrans(
    seasonal_cycle=seasonal_cycle, parameter_config=config
)

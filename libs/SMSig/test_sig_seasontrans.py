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

# Specify current directory create output directory if it does not exist
os.getcwd()

# %% [markdown]
# ## Read in observed data


# %%
def to_datetime(df, time_column, format="%Y-%m-%d %H:%M:%S"):
    df = df.copy()
    df[time_column] = pd.to_datetime(df[time_column], format=format)
    return df.set_index(time_column)


# %%

home_dir = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis"
data_dir = "data"
site = "Coweeta"

_data = pd.read_csv(
    r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\data\Coweeta\test_daily_2014_2018_sm_basinavg.csv"
)
_data = to_datetime(_data, "Time")
data = _data["Soil Moisture Content"]

data.head()

# %%
# Evaluate using seasonal soil moisture signature
sig_obs = SMSig(
    t=data.index.to_numpy(),
    sm=data.to_numpy(),
    plot_results=True,
    verbose=True,
)

# t_valley = sig_obs.calc_sinecurve()
# print(t_valley)

# %%
seasonal_cycle = pd.read_csv(
    os.path.join(home_dir, data_dir, site, "seasonal_cycle.csv"),
    parse_dates=["valley", "peak"],
)
seasonal_cycle


# %%

# Load the configuration
_parameter_config = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\data\Coweeta\seasonal_transition_config.json"
with open(_parameter_config, "r") as config_file:
    config = json.load(config_file)
season_trans_obs = sig_obs.calc_seasontrans(
    seasonal_cycle=seasonal_cycle, parameter_config=config
)

# %%

# Plot out the results
# df_obs = obs_synced
# df_sim = sim_synced
obs_label = "Observed"
sim_label = "Simulated"
obs_color = "#1f77b4"
sim_color = "#ff7f0e"
y_label = "Volumetric Soil Moisture Content [m^3/m^3]"
title = "Soil moisture and seasonal transition signatures"
fn = "timeseries.png"

# Relative values of SM
f2 = plt.figure(figsize=(30, 5))
# ax2 = f2.add_subplot(2,1,1)
# ax2.plot(sig_obs.tt.index, sig_obs.tt.values, alpha=1, label=obs_label, color=obs_color)
# ax2.plot(sig_sim.tt.index, sig_sim.tt.values, alpha=1, label=sim_label, color=sim_color)
# # ax2.plot(df["Time"].values, df_obs[var_name].values, alpha=1, label=obs_label, color=obs_color)
# # ax2.plot(df["Time"].values, df_sim[var_name].values, alpha=1, label=sim_label, color=sim_color)
# for i in range(len(start_dates_obs)):
#     ax2.axvline(x=start_dates_obs[i], color=obs_color, label=None, alpha=0.5, linestyle='-')
# for i in range(len(end_dates_obs)):
#     ax2.axvline(x=end_dates_obs[i], color=obs_color, label=None, alpha=0.5, linestyle='--')
# for i in range(len(start_dates_sim)):
#     ax2.axvline(x=start_dates_sim[i], color=sim_color, label=None, alpha=0.5, linestyle='-')
# for i in range(len(end_dates_sim)):
#     ax2.axvline(x=end_dates_sim[i], color=sim_color, label=None, alpha=0.5, linestyle='--')
# ax2.set_xlabel('Time')
# ax2.set_ylabel(y_label)
# ax2.set_title(title)
# ax2.legend()

ax3 = f2.add_subplot(2, 1, 2)
ax3.plot(
    sig_obs.tt.index,
    (sig_obs.tt.values - min(sig_obs.tt.values))
    / (max(sig_obs.tt.values) - min(sig_obs.tt.values)),
    alpha=1,
    label=obs_label,
    color=obs_color,
)
# ax3.plot(sig_sim.tt.index, (sig_sim.tt.values-min(sig_sim.tt.values))/(max(sig_sim.tt.values)-min(sig_sim.tt.values)), alpha=1, label=sim_label, color=sim_color)


def julian_to_datetime(jd):
    try:
        return datetime(1, 1, 1) + timedelta(days=jd - 1721425)
    except Exception as e:
        # Return np.nan if there's an error
        return np.nan


stard_dates_jd = np.concatenate([season_trans_obs[:, 0], season_trans_obs[:, 2]])
start_dates_obs = [julian_to_datetime(jd) for jd in stard_dates_jd]
end_dates_jd = np.concatenate([season_trans_obs[:1], season_trans_obs[:3]])
end_dates_obs = [julian_to_datetime(jd) for jd in stard_dates_jd]
for i in range(len(start_dates_obs)):
    ax3.axvline(
        x=start_dates_obs[i], color=obs_color, label=None, alpha=0.5, linestyle="-"
    )
for i in range(len(end_dates_obs)):
    ax3.axvline(
        x=end_dates_obs[i], color=obs_color, label=None, alpha=0.5, linestyle="--"
    )
# for i in range(len(start_dates_sim)):
#     ax3.axvline(x=start_dates_sim[i], color=sim_color, label=None, alpha=0.5, linestyle='-')
# for i in range(len(end_dates_sim)):
#     ax3.axvline(x=end_dates_sim[i], color=sim_color, label=None, alpha=0.5, linestyle='--')
ax3.set_xlabel("Time")
ax3.set_ylabel("Normalized soil moisture content")
# ax3.set_title(title)
ax3.legend()

# f2.savefig(os.path.join(out_path, fn), dpi=600)


# Save the results


# # %% [markdown]
# # ### Mahurangi

# # %%
# _data = pd.read_csv(
#     r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\data\Mahurangi\test_hourly_1998_2001_sm_basinavg.csv"
# )
# _data = to_datetime(_data, "Time", format=r"%m/%d/%Y %H:%M")
# data = _data["Soil Moisture Content"]

# data.head()

# # %%
# # Evaluate using seasonal soil moisture signature
# sig_obs = SMSig(
#     t=data.index.to_numpy(),
#     sm=data.to_numpy(),
#     plot_results=True,
#     verbose=True,
# )
# t_valley = sig_obs.calc_sinecurve()
# print(t_valley)

# # %%
# _t_valley_manual_input = pd.read_csv(
#     r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\data\Mahurangi\seasonal_cycel_valleys.csv",
#     header=None,
# )
# t_valley_manual_input = pd.to_datetime(_t_valley_manual_input[0])
# t_valley_manual_input

# # %%
# sig_obs = SMSig(
#     t=data.index.to_numpy(),
#     sm=data.to_numpy(),
#     plot_results=True,
#     verbose=True,
# )
# sig_obs.calc_seasontrans(t_valley=t_valley_manual_input)

# # %% [markdown]
# #

# %%

# %% Explanation from Readme document
	# - "Forcings" directory ----------------------------------------------
	# 	- NARR_day_9899.csv : Daily rainfall [m/s] on 1998-1999
	# 	- NARR_day_etp_9899.csv : Daily Potential Evapotranspiration [m/s] on 1998-1999
	# 	- NARR_day_9313.csv : Daily rainfall [m/s] on 1993-2013
	# 	- NARR_day_etp_9313.csv : Daily Potential Evapotranspiration [m/s] on 1993-2013

# PET data is daily 
# --> Extrapolate PET over the day by using temperature data (soil temp at 5cm is available), or flat rate (constant over the day)
# But if I use temp data, it's hard to define the min rate ... go with flat rate 
# How about mahurangi? -- hourly. intra-daily variation is already there
# Unit conversion
# --> m/s to m/hr, multiply by 60
# Time column is not clear ... seconds past from 1993-01-01? 
# https://zenodo.org/record/5289921

# %%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %%
filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\modelout_Picourlat_etal\Params\Forcings\NARR_day_etp_9313.csv"
df = pd.read_csv(filename, delim_whitespace=True, names=["time_sec", "PET_meter_per_sec"])

# %%
from datetime import datetime
df['datetime'] = pd.to_datetime(df['time_sec'], unit='s', origin='1993-01-01')
df.index = pd.date_range(periods=len(df), freq='86400s', start='1993-01-01')
df.tail()
# %%
# Convert PET [m/s] to [m/d], used in the CFE 
df['PET'] = df['PET_meter_per_sec'] * 60 * 60

# %%

# TODO: Resample to hourly timestep according to temperature data ()
df.set_index(df['datetime'])
df.index
# df_hourly = df.resample('H').ffill()

# %%
df_hourly = df.resample('H').ffill()
df_hourly.head()
# %%
# 0.00017 * 24 * 365 * 1000
# %%
filename = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp3\PET.csv"
df_hourly['PET'].to_csv(filename, index=True)
# %%
filename = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp3\PET.png"
# Plot and save the results
fig1, ax1 = plt.subplots()
df_hourly['PET'].plot(ax=ax1)
ax1.set_title(f'Little Washita')
ax1.set_xlabel("Time")
ax1.set_ylabel(r"[mm/hr]")
ax1.legend(loc='upper right')
fig1.savefig(filename)
# %%

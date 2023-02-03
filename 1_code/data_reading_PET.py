# %% Explanation from Readme document
	# - "Forcings" directory ----------------------------------------------
	# 	- NARR_day_9899.csv : Daily rainfall [m/s] on 1998-1999
	# 	- NARR_day_etp_9899.csv : Daily Potential Evapotranspiration [m/s] on 1998-1999
	# 	- NARR_day_9313.csv : Daily rainfall [m/s] on 1993-2013
	# 	- NARR_day_etp_9313.csv : Daily Potential Evapotranspiration [m/s] on 1993-2013

# PET data is daily 
# --> Extrapolate PET over the day by using temperature data (soil temp at 5cm is available), or flat rate (constant over the day)
# How about mahurangi? -- hourly. intra-daily variation is already there
# Unit conversion
# --> m/s to m/hr, multiply by 60
# Time column is not clear ... seconds past from 1993-01-01? 


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
df.tail()
# %%
# Convert PET [m/s] to [m/d], used in the CFE 
df['PET'] = df['PET_meter_per_sec'] * 60 * 60
df.head()
# %%

# TODO: Resample to hourly timestep according to temperature data 


# %%

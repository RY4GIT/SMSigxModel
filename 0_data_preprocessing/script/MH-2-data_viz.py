# %%
import pandas as pd
import numpy as np
import os

# %%
forcing_file_path = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\data\Mahurangi\mahurangi_1998_2001.csv"
data_file_path = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\data\Mahurangi\test_sm_basinavg.csv"

# %%
_forcing_data = pd.read_csv(forcing_file_path)
_forcing_data["time"] = pd.to_datetime(_forcing_data["time"])
_forcing_data.set_index("time", inplace=True)
forcing_data = _forcing_data
forcing_data.head()

# %%
_data = pd.read_csv(data_file_path)
_data["Time"] = pd.to_datetime(_data["Time"])
_data.set_index("Time", inplace=True)
data = _data
data.head()

# %%
sync_data = pd.merge_asof(data, forcing_data, left_index=True, right_index=True)
sync_data.head()

# %%
# Create a 3-by-1 plot using matplotlib
import matplotlib.pyplot as plt

fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

# Plot the approximated data in each subplot
sync_data["Rainfall"].plot(ax=axes[0])
axes[0].set_ylabel("Rainfall")

ax02 = axes[0].twinx()
sync_data["Flow"].plot(ax=ax02, color="orange")
ax02.set_ylabel("Flow")
ax02.invert_yaxis()

ax1 = axes[1]
sync_data["PET"].plot(ax=ax1)
ax1.set_ylabel("PET")

ax2 = axes[2]
sync_data["Soil Moisture Content"].plot(ax=ax2)
ax2.set_ylabel("VSWC [m3/m3]")

# Adjust layout to prevent overlap
plt.tight_layout()

plt.savefig(
    r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\data\Mahurangi\plot\input_data.png"
)

# Display the plot
plt.show()

# %%

# %%

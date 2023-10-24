# %%
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
input_path = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp3"
output_path = (
    r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\2_data_input\LittleWashita"
)
# %%
# df_alldata = pd.DataFrame()

if "df_alldata" in locals():
    del df_alldata

file_names = [
    "ars_VW05.csv",
    "ars_VW25.csv",
    "ars_VW45.csv",
    "ars_RAIN.csv",
    "PET.csv",
    "Q.csv",
]
variable_names = ["VW05", "VW25", "VW45", "RAIN", "PET", "Q"]


for i, filename in enumerate(file_names):
    variable_name = variable_names[i]

    print(filename)

    # # Read data
    df = pd.read_csv(os.path.join(input_path, filename))
    df["datetime"] = pd.to_datetime(df["datetime"])

    if variable_name != "Q":
        df["datetime"] = df["datetime"] + pd.Timedelta(hours=6)

    df.set_index(["datetime"], inplace=True)

    if variable_name != "Q":
        df.index = df.index.tz_localize("UTC")

    df.rename(columns={"0": variable_name}, inplace=True)

    # Sync them all
    if "df_alldata" in locals():
        df_alldata = df_alldata.join(df, how="outer")
    else:
        df_alldata = df

    del df

df_alldata.head()
# %%
df_alldata["Q"].plot()
# %%
df_alldata["VW05"].plot()
df_alldata["VW25"].plot()
df_alldata["VW45"].plot()
df_alldata["sm"] = df_alldata[["VW05", "VW25", "VW45"]].mean(axis=1).copy()

# %%
df_alldata["RAIN"] = df_alldata["RAIN"].copy() / 1000  # [mm/hr] --> [m/hr]
df_alldata["RAIN"].plot()
df_alldata["PET"].plot()

# %%
# Cut out
start_date = "2006-09-19"
end_date = "2012-09-19"
df_to_save = df_alldata[start_date:end_date].copy()
df_to_save = df_to_save.tz_convert(None)
# %% Plot

# %%
# Save
df_to_save_forcing = df_to_save[["RAIN", "PET"]].copy()
df_to_save_forcing.rename(columns={"RAIN": "precip_rate"}, inplace=True)
df_to_save_forcing.index.names = ["time"]
df_to_save_forcing.to_csv(
    os.path.join(output_path, "full", "forcing_hourly_2006_2012.csv")
)

# %%
area = 601 * 1e06
df_to_save_test = df_to_save[["RAIN", "Q", "sm"]].copy()
df_to_save_test.rename(
    columns={"sm": "Soil Moisture Content", "RAIN": "precip_rate", "Q": "Flow"},
    inplace=True,
)
df_to_save_test.index.names = ["Time"]
df_to_save_test["Time Step"] = df_to_save_test.reset_index().index
df_to_save_test["Direct Runoff"] = 0
df_to_save_test["GIUH Runoff"] = 0
df_to_save_test["Lateral Flow"] = 0
df_to_save_test["Base Flow"] = 0
df_to_save_test["Total Discharge"] = (
    df_to_save_test["Flow"] * area / 3600
)  # Is this correct? [m3/hr]?
df_to_save_test.to_csv(
    os.path.join(output_path, "full", "test_hourly_2006_2012_sm_basinavg.csv")
)

# %%

# Cut out short version
start_date = "2006-09-19"
end_date = "2007-09-19"
df_to_save_short_forcing = df_to_save_forcing[start_date:end_date].copy()
df_to_save_short_forcing.to_csv(
    os.path.join(output_path, "short", "forcing_hourly_2006_2012.csv")
)
df_to_save_short_test = df_to_save_test[start_date:end_date]
df_to_save_short_test.to_csv(
    os.path.join(output_path, "short", "test_hourly_2006_2012_sm_basinavg.csv")
)

# %%

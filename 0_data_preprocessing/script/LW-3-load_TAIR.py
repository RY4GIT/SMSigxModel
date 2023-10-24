# %%
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

# %%
# (0) Common: preparation
# Define the path to the folder containing the files
data_yrs = ["0609", "0911", "1113"]
temp_path = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp"
out_path = (
    r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\2_data_input\LittleWashita"
)

# %%
# (2) Then, concatinete sensors for the all period. Drop nan values, make time interval regular.

# Get a list of all files in the folder

temp_path_2 = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data_preprocessing\raw_data\LittleWashita\data_ars_temp2"
temp_path = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data_preprocessing\raw_data\LittleWashita\data_ars_temp"
# # Get the ones that has record from 2006
# sensors = [file.split('0609_a')[1][0:3] for file in file_list if (file != 'README') and ('0609' in file)]
n_yr = 0
file_list = os.listdir(temp_path)
print(file_list)
sensors = [
    file.split("_a")[1][0:3]
    for file in file_list
    if (file != "README") and ("0609" in file)
]
unique_sensors = list(set(sensors))
print(unique_sensors)
len(unique_sensors)

# %%
data_period_check = pd.DataFrame()
for i, sensor in enumerate(unique_sensors):
    print(f"Currently procesing sensor {i+1}/{len(unique_sensors)}")
    for n_yr in range(len(data_yrs)):
        file_path = os.path.join(
            temp_path, "ars" + data_yrs[n_yr] + "_a" + sensor + ".csv"
        )

        if os.path.exists(file_path):
            pass
        else:
            continue

        df = pd.read_csv(file_path, header=0)
        df["datetime"] = pd.to_datetime(df["datetime"])
        df = df.set_index(["datetime"])
        df.replace(-996, np.nan, inplace=True)
        df.replace(
            -998, np.nan, inplace=True
        )  # somehow the error flag in VW45 is not working
        if n_yr == 0:
            df_allyears = df
        else:
            df_allyears = pd.concat([df_allyears, df], axis=0)

    df_allyears.drop_duplicates()

    # Resample
    sampling_freq = "60min"
    df_allyears_daily = df_allyears.resample(sampling_freq).mean()

    file_name = "ars" + "_a" + str(sensor) + "_TAIR" + ".csv"
    file_path = os.path.join(temp_path_2, file_name)
    df_allyears_daily["TAIR"].to_csv(file_path, index=True)

    # Precipitation
    fig3, ax3 = plt.subplots()
    df_allyears_daily["TAIR"].plot(ax=ax3)
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Air Temp [C]")
    ax3.legend(loc="upper right")

    file_name = "ars" + "_a" + str(sensor) + "_TAIR.png"
    file_path = os.path.join(temp_path_2, file_name)
    fig3.savefig(os.path.join(file_path))

file_name = "ars_data_period_check.csv"
file_path = os.path.join(temp_path_2, file_name)
data_period_check.to_csv(file_path, index=False)

# %%
last_row = {
    "sensor_id": "shortest",
    "sensor_depth": "na",
    "start_date": data_period_check["start_date"].max(),
    "end_date": data_period_check["end_date"].min(),
}
data_period_check = data_period_check.append(last_row, ignore_index=True)
file_name = "ars_data_period_check_TAIR.csv"
file_path = os.path.join(temp_path_2, file_name)
data_period_check.to_csv(file_path, index=False)
# # %%
# data_period_check.tail()
# data_period_check.drop(data_period_check.tail(4).index, inplace = True)
# # %%
# data_period_check.tail()

"""
# %%
# %% (3) Put weights on it

thiessen_filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\gis\LittleWashita\Littel Washita\sm_sensors_thiessen.csv"
sitemeta_filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv"
temp_path_3 = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp3"


# %%
file_list = os.listdir(temp_path_2)
# print(file_list)
# sensors = [file.split('_a')[1][0:3] for file in file_list if (file != 'ars_data_period_check.csv')]
# unique_sensors = list(set(sensors))
# print(unique_sensors)
# len(unique_sensors)
# # data_variables = ['VW05', 'VW25', 'VW45', 'RAIN', 'TS05']
data_variables = ["RAIN"]
df_site_metadata_updated = pd.read_csv(
    r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv"
)
# %%

for j, data_variable in enumerate(data_variables):
    df_data = pd.DataFrame()
    df_area_fraction = pd.DataFrame()
    print(f"Currently processing {data_variable} data")

    for i, sensor in enumerate(unique_sensors):
        # Get the filename

        # target_filename = [filename for filename in file_list if (sensor in filename) and (data_variable in filename)]
        target_filename = "ars" + "_a" + str(sensor) + "_RAIN" + ".csv"
        print(target_filename)

        # Read the data
        df = pd.read_csv(os.path.join(temp_path_2, target_filename))
        df["datetime"] = pd.to_datetime(df["datetime"])
        df = df.set_index(["datetime"])
        df.rename({data_variable: sensor}, axis=1, inplace=True)
        df_data[sensor] = df

        # Get the area fraction
        target_area_fraction = df_site_metadata_updated[
            df_site_metadata_updated["stid"] == "A" + sensor
        ]["area_fraction_calculated_by_Ryoko"]
        df2 = pd.DataFrame([target_area_fraction.values] * len(df))
        df2.reindex_like(df)
        df2[np.isnan(df.values)] = np.nan
        df_area_fraction[sensor] = df2

    sum_fraction = df_area_fraction.sum(axis=1, skipna=True)
    sum_fraction.replace(0, np.nan, inplace=True)
    df_area_fraction_scaled = df_area_fraction.div(sum_fraction, axis="index")
    df_area_fraction_scaled

    df_data_weighted = df_data.values * df_area_fraction_scaled.values

    import matplotlib.pyplot as plt

    nparray_data_weighted_avg = np.nansum(df_data_weighted, axis=1)

    df_data_weighted_avg = pd.DataFrame(nparray_data_weighted_avg)
    df_data_weighted_avg.set_index(df_data.index, inplace=True)
    if data_variable != "RAIN":
        df_data_weighted_avg.replace(0, np.nan, inplace=True)

    file_name = "ars_" + data_variable + ".csv"
    file_path = os.path.join(temp_path_3, file_name)
    df_data_weighted_avg.to_csv(file_path, index=True)

    file_name = "ars_" + data_variable + ".png"
    file_path = os.path.join(temp_path_3, file_name)

    # Plot and save the results
    fig1, ax1 = plt.subplots()
    df_data_weighted_avg.plot(ax=ax1)
    ax1.set_title(f"Little Washita: sensor {data_variable}")
    ax1.set_xlabel("Time")
    ax1.set_ylabel(r"$\theta [m^3/m^3]$, [mm/hr]")
    ax1.legend(loc="upper right")
    fig1.savefig(os.path.join(file_path))

"""

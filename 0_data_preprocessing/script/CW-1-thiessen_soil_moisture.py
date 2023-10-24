# %%
# This code calculates the area fraction of Thiessen polygons and add the column to the metadata table
# For Little Washita, SM and rainfall sensors

# %%
import os
import pandas as pd
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt


def get_data(ncdf, variable_names_lst=[]):  # variable_names_lst=[] to get all the data
    # ws=watershed[watershed_number-1].split('-')[1]
    # print(ws)
    # File=ws+'_NetCDF.nc'
    # ncdf = nc4.Dataset(path +'/'+ File, 'r')
    keys_lst = list(ncdf.variables.keys())
    len_keys = len(keys_lst)
    print("Hydrometeorological variables in this watershed are:")
    all_variables = keys_lst[2 : int(len_keys / 2)]
    print(all_variables)
    if variable_names_lst == []:
        for i, var in enumerate(all_variables):
            print("Reading all variables")
            indexUnique = pd.date_range(
                str(pd.to_datetime(ncdf.variables["Datetime"][:][0])),
                str(pd.to_datetime(ncdf.variables["Datetime"][:][-1])),
            )
            var_df = pd.DataFrame(
                ncdf.variables[var][:],
                index=indexUnique,
                columns=[ncdf.variables[var].names],
            )
            if i == 0:
                former_df = var_df
            else:
                former_df = pd.concat([former_df, var_df], axis=1, join="outer")
    else:
        for i, var in enumerate(variable_names_lst):
            print(f"Reading {var}")
            indexUnique = pd.date_range(
                str(pd.to_datetime(ncdf.variables["Datetime"][:][0])),
                str(pd.to_datetime(ncdf.variables["Datetime"][:][-1])),
            )

            if ncdf.variables[var][:].shape[1] == 1:
                var_df = pd.DataFrame(
                    ncdf.variables[var][:],
                    index=indexUnique,
                    columns=[ncdf.variables[var].names],
                )
            else:
                var_df = pd.DataFrame(
                    ncdf.variables[var][:],
                    index=indexUnique,
                    columns=ncdf.variables[var].names,
                )

            if i == 0:
                former_df = var_df
            else:
                former_df = pd.concat([former_df, var_df], axis=1, join="outer")
    former_df.index.rename("DateTime", inplace=True)
    return former_df


# %%
thiessen_filename = (
    r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\gis\Coweeta\sm_thiessen_table.csv"
)
df_thiessen_area = pd.read_csv(thiessen_filename)
df_thiessen_area["area_fraction_calculated_by_Ryoko"] = (
    df_thiessen_area["Shape_Area"] / df_thiessen_area["Shape_Area"].sum()
)
print(df_thiessen_area.head())

# %%
sensors = df_thiessen_area.Sensor.values
unique_sensors = list(set(sensors))
print(unique_sensors)
len(unique_sensors)
# %%
# (3) Weight the rainfall and SM data according to Thiessen's polygon and take an avearge
temp_path = rf"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data_preprocessing\raw_data\Coweeta\temp"

# %%
# Read input data
data_dir = r"g:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data_preprocessing\raw_data\Coweeta\CHOSEN"
ncdf = nc4.Dataset(os.path.join(data_dir, f"Coweeta_NetCDF.nc"), "r")

df_sm = get_data(ncdf, ["SoilMoisture"])
ax = df_sm.filter(like="CWT102", axis=1).plot()
df_sm.filter(like="CWT202", axis=1).plot(ax=ax)
df_sm.filter(like="CWT302", axis=1).plot(ax=ax)
ax.set_xlim(["2015-06-01", "2018-12-31"])
ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

# %%

depths = ["30cm", "60cm"]
subsensors = ["A", "B"]
df_area_fraction = pd.DataFrame()
df_data_per_depth = pd.DataFrame()
for j, depth in enumerate(depths):
    for i, sensor in enumerate(unique_sensors):
        for k, subsensor in enumerate(subsensors):
            # Read the data correspponding to the sensor
            df_data = df_sm.filter(
                like=f"{sensor}_SoilMoisture_{depth}{subsensor}", axis=1
            ).copy()

            # Get the area fraction of the sensor of the interest
            target_area_fraction = df_thiessen_area[df_thiessen_area.Sensor == sensor][
                "area_fraction_calculated_by_Ryoko"
            ]

            # Create the array of area fraction, depending on the availability of the data
            df2 = pd.DataFrame([target_area_fraction.values] * len(df_data))
            df2.reindex_like(df_data)
            df2[np.isnan(df_data.values)] = np.nan
            df_area_fraction[f"{sensor}_{subsensor}"] = df2
            df_data_per_depth[f"{sensor}_{subsensor}"] = df_data

    # Get the sum of the fraction, and scale accordingly to get the true fraction
    sum_fraction = df_area_fraction.sum(axis=1, skipna=True)
    sum_fraction.replace(0, np.nan, inplace=True)
    df_area_fraction_scaled = df_area_fraction.div(sum_fraction, axis="index")
    df_area_fraction_scaled

    # Get the weighted average
    df_data_weighted = df_data_per_depth.values * df_area_fraction_scaled.values
    nparray_data_weighted_avg = np.nansum(df_data_weighted, axis=1)

    df_data_weighted_avg = pd.DataFrame(nparray_data_weighted_avg)
    df_data_weighted_avg.set_index(df_data_per_depth.index, inplace=True)
    df_data_weighted_avg.replace(0, np.nan, inplace=True)

    file_name = "sm_" + depth + "_avg.csv"
    file_path = os.path.join(temp_path, file_name)
    df_data_weighted_avg.to_csv(file_path, index=True)

    # %%
    # Plot and save the results

    fig1, ax1 = plt.subplots()
    df_data_weighted_avg.plot(ax=ax1)
    ax1.set_title(f"Coweeta: sensor {depth}")
    ax1.set_xlabel("Time")
    ax1.set_ylabel(r"$\theta [m^3/m^3]$, [mm/hr]")
    ax1.set_xlim("2015", "2019")
    ax1.legend(loc="upper right")

    file_name = "sm_" + depth + "_avg.png"
    file_path = os.path.join(temp_path, file_name)
    fig1.savefig(os.path.join(file_path))

    print(f"processed {depth} cm")


# %%

avg_data = pd.DataFrame()
for depth in depths:
    avg_data[depth] = pd.read_csv(
        os.path.join(temp_path, f"sm_{depth}_avg.csv"),
        index_col="DateTime",
        parse_dates=["DateTime"],
    )

avg_data["avg"] = avg_data.mean(axis=1)
avg_data["avg"].to_csv(os.path.join(temp_path, "sm_basin_avg.csv"), index=True)

# %%

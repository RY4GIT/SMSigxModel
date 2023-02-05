# %% 
# This code calculates the area fraction of Thiessen polygons and add the column to the metadata table 
# For Little Washita, SM and rainfall sensors

# %%
import os
import pandas as pd
import numpy as np
# %%
thiessen_filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\gis\LittleWashita\Littel Washita\sm_sensors_thiessen.csv"
sitemeta_filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv"

# # %%
# df_thiessen_area = pd.read_csv(thiessen_filename)
# print(df_thiessen_area.head())
# df_site_metadata = pd.read_csv(sitemeta_filename)
# print(df_site_metadata.head())
# # %%
# df_site_metadata_updated = pd.merge(df_site_metadata, df_thiessen_area[['stid', 'Shape_Area']], how='left', on=['stid'])
# df_site_metadata_updated.head()
# # %%
# df_site_metadata_updated['area_fraction_calculated_by_Ryoko'] = df_site_metadata_updated['Shape_Area']/df_site_metadata_updated['Shape_Area'].sum()
# # %%
# df_site_metadata_updated['area_fraction_calculated_by_Ryoko'].sum()
# # %%
# df_site_metadata_updated.to_csv(sitemeta_filename)
# # %%

# %%
# (3) Weight the rainfall and SM data according to Thiessen's polygon and take an avearge
temp_path_2 = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp2'
temp_path_3 = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp3'

# %%
file_list = os.listdir(temp_path_2)
print(file_list)
sensors = [file.split('_a')[1][0:3] for file in file_list if (file != 'ars_data_period_check.csv')]
unique_sensors = list(set(sensors))
print(unique_sensors)
len(unique_sensors)
data_variables = ['VW05', 'VW25', 'VW45', 'RAIN', 'TS05']
df_site_metadata_updated = pd.read_csv(r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv")
# %%

for j, data_variable in enumerate(data_variables):
    df_data = pd.DataFrame()
    df_area_fraction = pd.DataFrame()
    print(f'Currently processing {data_variable} data')
    
    for i, sensor in enumerate(unique_sensors):
        # Get the filename 
        target_filename = [filename for filename in file_list if (sensor in filename) and (data_variable in filename)]
        print(target_filename)

        # Read the data 
        df = pd.read_csv(os.path.join(temp_path_2, target_filename[0]))
        df['datetime'] = pd.to_datetime(df['datetime'])
        df = df.set_index(['datetime'])     
        df.rename({data_variable: sensor}, axis=1, inplace=True)
        df_data[sensor] = df
        
        # Get the area fraction
        target_area_fraction = df_site_metadata_updated[df_site_metadata_updated['stid'] == 'A'+ sensor]['area_fraction_calculated_by_Ryoko']
        df2 = pd.DataFrame([target_area_fraction.values] * len(df))
        df2.reindex_like(df)
        df2[np.isnan(df.values)] = np.nan
        df_area_fraction[sensor] = df2

    sum_fraction = df_area_fraction.sum(axis=1, skipna=True)
    sum_fraction.replace(0, np.nan, inplace=True)
    df_area_fraction_scaled = df_area_fraction.div(sum_fraction, axis='index')
    df_area_fraction_scaled


    df_data_weighted = df_data.values * df_area_fraction_scaled.values

    import matplotlib.pyplot as plt
    nparray_data_weighted_avg = np.nansum(df_data_weighted, axis=1)


    df_data_weighted_avg = pd.DataFrame(nparray_data_weighted_avg) 
    df_data_weighted_avg.set_index(df_data.index, inplace=True)
    if data_variable != 'RAIN':
        df_data_weighted_avg.replace(0, np.nan, inplace=True)

    file_name = 'ars_' + data_variable + '.csv'
    file_path = os.path.join(temp_path_3, file_name)
    df_data_weighted_avg.to_csv(file_path, index=True)

    file_name = 'ars_' + data_variable + '.png'
    file_path = os.path.join(temp_path_3, file_name)

    # Plot and save the results
    fig1, ax1 = plt.subplots()
    df_data_weighted_avg.plot(ax=ax1)
    ax1.set_title(f'Little Washita: sensor {data_variable}')
    ax1.set_xlabel("Time")
    ax1.set_ylabel(r"$\theta [m^3/m^3]$, [mm/hr]")
    ax1.legend(loc='upper right')
    fig1.savefig(os.path.join(file_path))


# %%

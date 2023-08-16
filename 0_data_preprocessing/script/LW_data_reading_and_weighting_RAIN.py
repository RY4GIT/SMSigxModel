# %% [markdown]
# This code read and calculate weights on the rainfall dataset. I realized the errors in my previous code later, so the previous code may not work now...

# %%
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# %%
data_yrs = ['0609', '0911', '1113']
temp_path = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp'
out_path = r'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\2_data_input\LittleWashita'


# %% [markdown]
# ## (1) First, concatinate all data within a folder and store them in the temporary folder

# %%
for n_yr in range(len(data_yrs)):
    # Extract the dates and sensor numbers from the file names
    folder_path = rf'G:|\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_{data_yrs[n_yr]}'
    # Get a list of all files in the folder
    file_list = os.listdir(folder_path)
    # file_list[1].split('a')[1][0:3]

    dates = [file.split('a')[0] for file in file_list if file != 'README']
    sensors = [file.split('a')[1][0:3] for file in file_list if file != 'README']
    
    unique_dates = list(set(dates))
    unique_sensors = list(set(sensors))
    
    # Create lists to store the dataframes for each sensor
    sensor_dataframes = {}

    # for sensor in unique_sensors:
    #     # Create an empty dataframe for the sensor
    #     sensor_dataframes[sensor] = pd.DataFrame()
    for i, sensor in enumerate(unique_sensors):
        sensor_dataframes[sensor] = pd.DataFrame()
        print(f'Currently processing sensor {i+1}/{len(unique_sensors)}')
        
        for date in tqdm(unique_dates):
            # Extract the date and sensor number from the file name
            file_name = os.path.join(folder_path, date + 'a' + sensor + '.txt')
            
            # Read the file as a pandas dataframe
            file_path = os.path.join(folder_path, file_name)
            
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace=True, header=3)
                
                df['date'] = pd.to_datetime(date, format='%Y%m%d')
                
                mask = df['HR'] == 24
                df.loc[mask, 'date'] = df.loc[mask, 'date'] + pd.Timedelta(days=1)
                df.loc[mask, 'HR'] = '00'
                
                df['datetime'] = pd.to_datetime(df['date'].astype(str) + ' ' + df['HR'].astype(str).str.zfill(2) + ':'+ df['MN'].astype(str).str.zfill(2), format='%Y-%m-%d %H:%M')
                
                df.loc[df["QRN"] != "g", "RAIN"] = np.nan
                rain_at_the_start = df["RAIN"][0].copy()
                df["RAIN"] = df["RAIN"].diff() 
                # Convert cumulative to original
                df["RAIN"][0] = rain_at_the_start

                df = df[["datetime", "RAIN"]]

                # Concatenate the dataframe into the dataframe for the sensor
                sensor_dataframes[sensor] = pd.concat([sensor_dataframes[sensor], df], axis=0, ignore_index=True)
            else:
                pass
        
        # Save the dataframe for the sensor to a CSV file
        file_name = 'ars' + data_yrs[n_yr] + '_a' + sensor + '_rain.csv'
        file_path = os.path.join(temp_path, file_name)
        sensor_dataframes[sensor].to_csv(file_path, index=False)

# %% [markdown]
# # (2) Then, concatinete sensors for the all period. Drop nan values, make time interval regular. 

# %%

temp_path_2 = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp2'
temp_path = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp'
# # Get the ones that has record from 2006
# sensors = [file.split('0609_a')[1][0:3] for file in file_list if (file != 'README') and ('0609' in file)]
n_yr = 0
file_list = os.listdir(temp_path)
print(file_list)
sensors = [file.split('_a')[1][0:3] for file in file_list if (file != 'README') and ('0609' in file)]
unique_sensors = list(set(sensors))
print(unique_sensors)
len(unique_sensors)


# %%
data_period_check = pd.DataFrame()
for i, sensor in enumerate(unique_sensors):
    print(f'Currently procesing sensor {i+1}/{len(unique_sensors)}')
    for n_yr in range(len(data_yrs)):
        file_path = os.path.join(temp_path, 'ars' + data_yrs[n_yr] + '_a' + sensor + '_rain.csv')
        
        if os.path.exists(file_path):
            pass
        else:
            continue
        
        df = pd.read_csv(file_path, header=0)
        df['datetime'] = pd.to_datetime(df['datetime'])
        df = df.set_index(['datetime'])
        df.replace(-996, np.nan, inplace=True)
        df.replace(-998, np.nan, inplace=True)# somehow the error flag in VW45 is not working 
        if n_yr == 0:
            df_allyears = df
        else: 
            df_allyears = pd.concat([df_allyears, df], axis=0)
        
    df_allyears.drop_duplicates() 

    # Resample
    sampling_freq = '60min'
    df_allyears = df_allyears.fillna(0)
    df_allyears_hourly = df_allyears.resample(sampling_freq).sum() # Conversion
    
    file_name = 'ars' + '_a' + str(sensor) + '_RAIN' + '.csv'
    file_path = os.path.join(temp_path_2, file_name)
    df_allyears_hourly['RAIN'].to_csv(file_path, index=True)
    
    # Precipitation
    fig2, ax2 = plt.subplots()
    df_allyears_hourly['RAIN'].plot(ax=ax2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Precip [mm/hr]")
    ax2.legend(loc='upper right')
    
    file_name = 'ars' + '_a' + str(sensor) + '_RAIN.png'
    file_path = os.path.join(temp_path_2, file_name)
    fig2.savefig(os.path.join(file_path))


# %% (3) Put weights on it

thiessen_filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\gis\LittleWashita\Littel Washita\sm_sensors_thiessen.csv"
sitemeta_filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv"
temp_path_3 = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp3'


# %%
file_list = os.listdir(temp_path_2)
# print(file_list)
# sensors = [file.split('_a')[1][0:3] for file in file_list if (file != 'ars_data_period_check.csv')]
# unique_sensors = list(set(sensors))
# print(unique_sensors)
# len(unique_sensors)
# # data_variables = ['VW05', 'VW25', 'VW45', 'RAIN', 'TS05']
data_variables = ['RAIN']
df_site_metadata_updated = pd.read_csv(r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv")
# %%

for j, data_variable in enumerate(data_variables):
    df_data = pd.DataFrame()
    df_area_fraction = pd.DataFrame()
    print(f'Currently processing {data_variable} data')
    
    for i, sensor in enumerate(unique_sensors):
        # Get the filename 
    
        # target_filename = [filename for filename in file_list if (sensor in filename) and (data_variable in filename)]
        target_filename = 'ars' + '_a' + str(sensor) + '_RAIN' + '.csv'
        print(target_filename)

        # Read the data 
        df = pd.read_csv(os.path.join(temp_path_2, target_filename))
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

# %%
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

# %%
# (0) Common: preparation
# Define the path to the folder containing the files
data_yrs = ['0609', '0911', '1113']
temp_path = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp'
out_path = r'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\2_data_input\LittleWashita'

# %%
# (1) First, concatinate all data within a folder and store them in the temporary folder
for n_yr in range(len(data_yrs)):
    # Extract the dates and sensor numbers from the file names
    folder_path = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_{data_yrs[n_yr]}'
    # Get a list of all files in the folder
    file_list = os.listdir(folder_path)
    # file_list[1].split('a')[1][0:3]

    dates = [file.split('a')[0] for file in file_list if file != 'README']
    sensors = [file.split('a')[1][0:3] for file in file_list if file != 'README']
    
    # %%
    unique_dates = list(set(dates))
    unique_sensors = list(set(sensors))
    # %%
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
                df.loc[df["QTS5"] != "g", "TS05"] = np.nan
                df.loc[df["QVW5"] != "g", "VW05"] = np.nan
                df.loc[df["QVW25"] != "g", "VW25"] = np.nan
                df.loc[df["QVW45"] != "g", "VW45"] = np.nan
                
                if 'TAIR' not in df.columns:
                    # Add "TAIR" column filled with NaN if it doesn't exist
                    df = df.assign(TAIR=np.full(len(df), np.nan))
                else:
                    df.loc[df["TAIR"] != "g", "QTA"] = np.nan

                df = df[["datetime", "RAIN", "TAIR", "QTS5", "VW05", "VW25", "VW45"]]

                # Concatenate the dataframe into the dataframe for the sensor
                sensor_dataframes[sensor] = pd.concat([sensor_dataframes[sensor], df], axis=0, ignore_index=True)
            else:
                pass
        
        # Save the dataframe for the sensor to a CSV file
        file_name = 'ars' + data_yrs[n_yr] + '_a' + sensor + '.csv'
        file_path = os.path.join(temp_path, file_name)
        sensor_dataframes[sensor].to_csv(file_path, index=False)
    
# %%
# (2) Then, concatinete sensors for the all period. Drop nan values, make time interval regular. 

# Get a list of all files in the folder

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
        file_path = os.path.join(temp_path, 'ars' + data_yrs[n_yr] + '_a' + sensor + '.csv')
        
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
    
    # if os.path.exists(file_path):
    #     pass
    # else:
    #     continue
        
    df_allyears.drop_duplicates() 

    # Resample
    sampling_freq = '60min'
    df_allyears_daily = df_allyears.resample(sampling_freq).mean()

    # Data check
    for sensor_depth in ['VW05', 'VW25', 'VW45']:
        new_row = {'sensor_id':sensor, 'sensor_depth': sensor_depth, 'start_date': df_allyears_daily[sensor_depth].first_valid_index(), 'end_date': df_allyears_daily[sensor_depth].last_valid_index()}
        data_period_check = data_period_check.append(new_row, ignore_index=True)

    # Save the results
    for sensor_depth in ['VW05', 'VW25', 'VW45']:
        file_name = 'ars' + '_a' + str(sensor) + '_d' + sensor_depth + '.csv'
        file_path = os.path.join(temp_path_2, file_name)
        df_allyears_daily[sensor_depth].to_csv(file_path, index=True)
        
    file_name = 'ars' + '_a' + str(sensor) + '_RAIN' + '.csv'
    file_path = os.path.join(temp_path_2, file_name)
    df_allyears_daily['RAIN'].to_csv(file_path, index=True)
    
    file_name = 'ars' + '_a' + str(sensor) + '_TAIR' + '.csv'
    file_path = os.path.join(temp_path_2, file_name)
    df_allyears_daily['RAIN'].to_csv(file_path, index=True)
    
    file_name = 'ars' + '_a' + str(sensor) + '_TS05' + '.csv'
    file_path = os.path.join(temp_path_2, file_name)
    df_allyears_daily['TS05'].to_csv(file_path, index=True)
    
    # Plot and save the results
    fig1, ax1 = plt.subplots()
    df_allyears_daily[['VW05', 'VW25', 'VW45']].plot(ax=ax1)
    ax1.set_title(f'Little Washita: sensor a{sensor}')
    ax1.set_xlabel("Time")
    ax1.set_ylabel(r"$\theta [m^3/m^3]$")
    ax1.legend(loc='upper right')
    
    file_name = 'ars' + '_a' + str(sensor) + '_SM.png'
    file_path = os.path.join(temp_path_2, file_name)
    fig1.savefig(os.path.join(file_path))

    # Precipitation
    fig2, ax2 = plt.subplots()
    df_allyears_daily['RAIN'].plot(ax=ax2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Precip [mm/hr]")
    ax2.legend(loc='upper right')
    
    file_name = 'ars' + '_a' + str(sensor) + '_RAIN.png'
    file_path = os.path.join(temp_path_2, file_name)
    fig2.savefig(os.path.join(file_path))
    
    # Precipitation
    fig3, ax3 = plt.subplots()
    df_allyears_daily['TS05'].plot(ax=ax3)
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Air Temp [C]")
    ax3.legend(loc='upper right')
    
    file_name = 'ars' + '_a' + str(sensor) + '_TS05.png'
    file_path = os.path.join(temp_path_2, file_name)
    fig3.savefig(os.path.join(file_path))

file_name = 'ars_data_period_check.csv'
file_path = os.path.join(temp_path_2, file_name)
data_period_check.to_csv(file_path, index=False)

# %%
last_row = {'sensor_id':'shortest', 'sensor_depth': 'na', 'start_date': data_period_check['start_date'].max(), 'end_date': data_period_check['end_date'].min()}
data_period_check = data_period_check.append(last_row, ignore_index=True)
file_name = 'ars_data_period_check.csv'
file_path = os.path.join(temp_path_2, file_name)
data_period_check.to_csv(file_path, index=False)

# # %%
# data_period_check.tail()
# data_period_check.drop(data_period_check.tail(4).index, inplace = True)
# # %%
# data_period_check.tail()





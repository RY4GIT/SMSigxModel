# %%
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
# %%

# Define the path to the folder containing the files
data_yrs = ['0609', '0911', '1113']
n_yr = 2

folder_path = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_{data_yrs[n_yr]}'
out_path = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\data_ars_temp'
# Get a list of all files in the folder
file_list = os.listdir(folder_path)
# file_list[1].split('a')[1][0:3]
# %%
# Extract the dates and sensor numbers from the file names
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
        file_path = os.path.join(folder_path, date + 'a' + sensor + '.txt')
        
        # Read the file as a pandas dataframe
        file_path = os.path.join(folder_path, file_path)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, delim_whitespace=True, header=3)
            
            df['date'] = pd.to_datetime(date, format='%Y%m%d')
            
            mask = df['HR'] == 24
            df.loc[mask, 'date'] = df.loc[mask, 'date'] + pd.Timedelta(days=1)
            df.loc[mask, 'HR'] = '00'
            
            df['datetime'] = pd.to_datetime(df['date'].astype(str) + ' ' + df['HR'].astype(str).str.zfill(2) + ':'+ df['MN'].astype(str).str.zfill(2), format='%Y-%m-%d %H:%M')
            
            df.loc[df["QRN"] != "g", "RAIN"] = np.nan
            df.loc[df["QVW5"] != "g", "VW05"] = np.nan
            df.loc[df["QVW25"] != "g", "VW25"] = np.nan
            df.loc[df["QVW45"] != "g", "VW45"] = np.nan
            
            df = df[["datetime", "RAIN", "VW05", "VW25", "VW45"]]

            # Concatenate the dataframe into the dataframe for the sensor
            sensor_dataframes[sensor] = pd.concat([sensor_dataframes[sensor], df], axis=0, ignore_index=True)
        else:
            pass
    
    # Save the dataframe for the sensor to a CSV file
    file_name = 'ars' + data_yrs[n_yr] + '_a' + sensor + '.csv'
    file_path = os.path.join(out_path, file_name)
    sensor_dataframes[sensor].to_csv(file_path, index=False)
    
# %%

# %%

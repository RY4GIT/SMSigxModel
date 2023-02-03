# %% 
# This code calculates the area fraction of Thiessen polygons and add the column to the metadata table 
# For Little Washita, SM and rainfall sensors

# %%
import os
import pandas as pd
import numpy as np
# %%
filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\gis\LittleWashita\Littel Washita\sm_sensors_thiessen.csv"
df_thiessen_area = pd.read_csv(filename)
print(df_thiessen_area.head())

filename = r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv"
df_site_metadata = pd.read_csv(filename)
print(df_site_metadata.head())
# %%
df_site_metadata_updated = pd.merge(df_site_metadata, df_thiessen_area[['stid', 'Shape_Area']], how='left', on=['stid'])
df_site_metadata_updated.head()
# %%

df_site_metadata_updated['area_fraction_calculated_by_Ryoko'] = df_site_metadata_updated['Shape_Area']/df_site_metadata_updated['Shape_Area'].sum()
# %%
df_site_metadata_updated['area_fraction_calculated_by_Ryoko'].sum()
# %%
df_site_metadata_updated.to_csv(r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita\site_metadata_ars.csv")
# %%

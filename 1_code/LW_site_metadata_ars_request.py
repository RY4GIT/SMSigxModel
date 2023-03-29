# %%
# http://api.mesonet.org/index.php/meta/ars_sites/network/all/status/all
import requests

url = "http://api.mesonet.org/index.php/meta/ars_sites/network/all/status/all"

response = requests.get(url)

if response.status_code == 200:
    data = response.json()
    # Do something with the data
else:
    print("Request failed with status code:", response.status_code)

# dict_keys(['response', 'api_version', 'params'])


# %%
data['response']

# %%
data['params']
# %%
import pandas as pd
site_metadata = pd.DataFrame()
for i in range(len(data['response'])):
    df = pd.DataFrame([data['response'][i]])
    site_metadata = pd.concat([site_metadata, df], axis=0, ignore_index=True)
# %%
site_metadata.head()
site_metadata.sort_values(by='status')

import os
out_path = rf'G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\0_data\Little Washita'
file_name = 'site_metadata_ars.csv'
file_path = os.path.join(out_path, file_name)
site_metadata.to_csv(file_path, index=False)

# %%

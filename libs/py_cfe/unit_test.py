# %% [markdown]
# ### This jupyter notebook is to unit-test CFE model run

# %% [markdown]
# ## Preparation

# %%
# Import modules
import sys
import json
from bmi_cfe import BMI_CFE

# %%
# Input: Mahurangi data for 3 yrs
# forcing: Mahurangi/mahurangi_1998_2001.csv
# observed: test_sm_basinavg

# Input: Little Washita data for 12 yrs
# forcing: LittleWashita/little_washita_2006_2012.csv
# observed: test_sm_basinavg


# Parameter: only limited number of parameters

input_json = "./data/unit_test/config_cfe.json"
params = {
    "forcing_file": "G:/Shared drives/Ryoko and Hilary/SMSigxModel/analysis/data/LittleWashita/little_washita_2006_2012.csv",
    "catchment_area_km2": 601,
    "soil_params": {
        "bb": 6.591999812593777,
        "satdk": 0.00001,
        "slop": 0.0092968319200453,
        "satpsi": 0.1736243088137391,
        "smcmax": 0.5,
        "wltsmc": 0.0129467948981879,
        "D": 2,
    },
    "alpha_fc": 0.3083778619066096,
    "max_gw_storage": 0.1061405978955891,
    "Cgw": 0.0009875396369054,
    "expon": 3.153967939163099,
    "K_lf": 0.5731958888031214,
    "refkdt": 10000,
    "K_nash": 0.6310887218737914,
    "trigger_z_fact": 0.6381182019616936,
    "num_nash_storage": 2,
    "nash_storage": [0.0, 0.0],
    "giuh_ordinates": [0.1, 0.4, 0.2, 0.2, 0.1],
    "stand_alone": 1,
    "unit_test": 1,
    "compare_results_file": "G:/Shared drives/Ryoko and Hilary/SMSigxModel/analysis/data/LittleWashita/test_sm_basinavg.csv",
}

with open(input_json, "w") as outfile:
    json.dump(params, outfile, indent=4)

with open(input_json) as outfile:
    loaded_data = json.load(outfile)

# %% [markdown]
# ## Run the model

# %%
cfe_instance = BMI_CFE(cfg_file=input_json, verbose=True)
cfe_instance.initialize()
cfe_instance.run_unit_test(plot=True, print_fluxes=False, warm_up=True)
cfe_instance.finalize(print_mass_balance=True)

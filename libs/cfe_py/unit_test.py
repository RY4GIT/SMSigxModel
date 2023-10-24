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
    "forcing_file": "G:/Shared drives/Ryoko and Hilary/SMSigxModel/analysis/data/Coweeta/forcing_daily_2014_2018.csv",
    "catchment_area_km2": 0.1210,
    "soil_params": {
        "bb": 4.74,
        "satdk": 2.45e-06,
        "slop": 0.147,
        "satpsi": 0.263,
        "smcmax": 0.5,
        "wltsmc": 0.35,
        "D": 2,
    },
    "max_gw_storage": 50,
    "Cgw": 1,
    "expon": 1.75,
    "K_lf": 0.005,
    "K_nash": 0.1,
    "num_nash_storage": 2,
    "giuh_ordinates": [1.0],
    "trigger_z_fact": 0.5,
    "alpha_fc": 0.18,
    "refkdt": 3,
    "stand_alone": 1,
    "unit_test": 1,
    "time_step_size": 86400,
    "compare_results_file": "G:/Shared drives/Ryoko and Hilary/SMSigxModel/analysis/data/Coweeta/test_daily_2014_2018_sm_basinavg.csv",
}

with open(input_json, "w") as outfile:
    json.dump(params, outfile, indent=4)

with open(input_json) as outfile:
    loaded_data = json.load(outfile)

# %% [markdown]
# ## Run the model

# %%
cfe_instance = BMI_CFE(cfg_file=input_json, verbose=False)
cfe_instance.initialize()
cfe_instance.run_unit_test(
    plot=True,
    print_fluxes=True,
    plot_lims=list(range(1, 1034)),
    warm_up=True,
    warmup_offset=350,
    warmup_iteration=10,
)
cfe_instance.finalize(print_mass_balance=True)

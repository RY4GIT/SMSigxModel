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
        "bb": 4.74,
        "satdk": 1e-02,
        "slop": 1,
        "satpsi": 0.141,
        "smcmax": 0.4,
        "wltsmc": 0.05,
        "D": 2,
        "exponent_secondary": 0.71689,
    },
    "allow_percolation_below_threshold": 0,
    "alpha_fc": 0.18,
    "max_gw_storage": 1,
    "Cgw": 1e-04,
    "expon": 0.60668 * 2,
    "gw_scheme": "Linear",
    "K_lf": 1e-03,
    "refkdt": 500,
    "K_nash": 0.01,
    "trigger_z_fact": 0.25,
    "num_nash_storage": 250,
    "giuh_ordinates": [0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.4, 0.2, 0.2, 0.1],
    "stand_alone": 1,
    "unit_test": 1,
    # "revap_factor": 0.05,
    "compare_results_file": "G:/Shared drives/Ryoko and Hilary/SMSigxModel/analysis/data/LittleWashita/test_sm_basinavg.csv",
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
    print_fluxes=False,
    plot_lims=list(range(1, 52600)),
    warm_up=True,
    warmup_offset=12000,
    warmup_iteration=10,
)
cfe_instance.finalize(print_mass_balance=True)

import spotpy
import pandas as pd
import json


def salib_cfe(X, param_names, myCFE):

    # write the randomly-generated parameters to the config json file
    with open(myCFE.cfg_file) as data_file:
        cfe_cfg = json.load(data_file)

    for i, param_name in enumerate(param_names):
        if param_name in ['depth', 'bb', 'mult', 'satdk', 'satpsi', 'slop', 'smcmax', 'wltsmc', 'D']:
            cfe_cfg["soil_params"][param_name] = X[i]
        else:
            cfe_cfg[param_name] = X[i]

    with open(myCFE.cfg_file, 'w') as out_file:
        json.dump(cfe_cfg, out_file)

    # Here the model is actualy started with a unique parameter combination that it gets from spotpy for each time the model is called
    myCFE.initialize()
    myCFE.run_unit_test(plot=False, print_fluxes=False)

    sim = myCFE.cfe_output_data[["Time", "Total Discharge"]]
    sim["Time"] = pd.to_datetime(sim["Time"])
    sim = sim.set_index("Time")
    # return sim

    # Get the comparison data
    data = myCFE.unit_test_data
    obs = data[["Time", "Total Discharge"]]
    obs["Time"] = pd.to_datetime(obs["Time"])
    obs = obs.set_index("Time")
    # return obs

    df = pd.merge_asof(sim, obs, on = "Time")
    sim_synced = df["Total Discharge_x"]
    obs_synced = df["Total Discharge_y"]

    # Calculate objective metrics
    like = spotpy.objectivefunctions.nashsutcliffe(obs_synced, sim_synced)

    return like




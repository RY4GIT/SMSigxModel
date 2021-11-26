import spotpy
import pandas as pd
import json

def salib_cfe(X, param_names, myCFE):

    # write the randomly-generated parameters to the config json file
    with open(myCFE.cfg_file) as data_file:
        cfe_cfg = json.load(data_file)

    for i, param_name in enumerate(param_names):
        cfe_cfg["soil_params"][param_name] = X[i]

    with open(myCFE.cfg_file, 'w') as out_file:
        json.dump(cfe_cfg, out_file)

    # Here the model is actualy started with a unique parameter combination that it gets from spotpy for each time the model is called
    data = myCFE.run_unit_test(plot_results=False)
    sim = data[["Time", "Total Discharge"]]
    sim["Time"] = pd.to_datetime(sim["Time"])
    return sim

    # Get the comparison data
    data = myCFE.load_unit_test_data()
    obs = data[["Time", "Total Discharge"]]
    obs["Time"] = pd.to_datetime(obs["Time"])
    return obs

    df = pd.merge_asof(simulation, evaluation, on = "Time")
    sim_synced = df["Total Discharge_x"]
    obs_synced = df["Total Discharge_y"]

    # Calculate objective metrics
    like = spotpy.objectivefunctions.nashsutcliffe(obs_synced, sim_synced)

    return like




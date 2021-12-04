import spotpy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.sample import morris as morris_s
from SALib.analyze import morris as morris_a


class SALib_CFE():

    def __init__(self, cfe_instance=None, problem=None, SAmethod=None):
        self.cfe_instance = cfe_instance
        self.problem = problem
        self.SAmethod = SAmethod

    def run(self):
        if self.SAmethod == "Sobol":
            # sample
            n = 10
            param_values = saltelli.sample(self.problem, n, calc_second_order=True)

            # run a model
            Y = run_cfes(
                problem = self.problem,
                cfe_instance = self.cfe_instance,
                param_values=param_values,
                nrun = n*(2*self.problem['num_vars']+2)
            )

            Si = sobol.analyze(self.problem, param_values, Y, calc_second_order=True, print_to_console=False)
            # TypeError: analyze() got multiple values for argument 'calc_second_order'
            print(Si)

        if self.SAmethod == "Morris":
            # sample
            iteration = 10
            n_levels = 4
            param_values = morris_s.sample(self.problem, iteration, num_levels=n_levels)

            # run a model
            Y = run_cfes(
                problem = self.problem,
                cfe_instance = self.cfe_instance,
                param_values=param_values,
                nrun = iteration * n_levels
            )

            # evaluation
            Si = morris_a.analyze(self.problem, param_values, Y, print_to_console=True)
            print(Si['mu'])

    def plot(self):
        None


def run_cfes(problem, cfe_instance, param_values, nrun):
    Y = np.zeros([param_values.shape[0]])
    for i, X in enumerate(param_values):
        print('{} of {}'.format(i, nrun))
        Y[i] = salib_cfe_interface(X, problem['names'], cfe_instance)
    return Y


def salib_cfe_interface(X, param_names, myCFE):

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




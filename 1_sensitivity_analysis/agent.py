
import os
import warnings
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from SALib.sample import morris as morris_s
from SALib.analyze import morris as morris_a

from tqdm import tqdm
from model import CFEmodel


def read_SALib_config(config_SALib_path):
    """Reads configuration setting for SALib sensitivity analysis from text file."""

    # Read the data
    df = pd.read_csv(config_SALib_path)
    df_param_to_calibrate = df[df['calibrate']==1]

    # Convert it to dictionary
    config_SALib = {
        'num_vars': df_param_to_calibrate.shape[0],
        'names': df_param_to_calibrate['name'].tolist(),
        'bounds': df_param_to_calibrate[['lower_bound', 'upper_bound']].values.tolist()
    }

    for i, bounds in enumerate(config_SALib['bounds']):
        if bounds[1] <= bounds[0]:
            warnings.warn(f"The upper bound is smaller than the lower bound for the parameter {config_SALib['names'][i]}")
    return config_SALib
        
class Agent_SALib_CFE():

    def __init__(self, config=None):
        
        self.config = config
        self.problem = read_SALib_config(config['PATHS']['salib_config'])
        
        self.out_path = os.path.join(config['PATHS']['homedir'], 'results', self.config['DATA']['site'])
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)
            
    def run(self):
        
        runtype = self.config['SALib']['method']
        
        if runtype == "Morris":
            self.run_Morris()
        else:
            print(f"Invalid runtype: {runtype}")
    
    def run_Morris(self):
        
        # Sample parameters 
        N = int(self.config['Morris']['N'])
        n_levels = int(self.config['Morris']['n_levels'])
        self.sampled_params = morris_s.sample(self.problem, N=N, num_levels=n_levels)
        
        # Define number of runs 
        nrun = N * (self.problem['num_vars']+1)
        print(f'Total runs: {nrun} \n Number of analyzed parameters: {self.problem["num_vars"]}\n')
        
        # Initialize output array 
        self.Y = np.zeros([self.sampled_params.shape[0]])
        
        # Run the simulations and evaluation
        for i, X in tqdm(enumerate(self.sampled_params)):
            self.model = CFEmodel(config=self.config, problem=self.problem, X=X)
            self.model.run()
            self.Y[i] = self.model.evaluate()

        # Runo morris
        print("### Results ###")
        self.Si = morris_a.analyze(self.problem, self.sampled_params, self.Y, print_to_console=True)

        # Output the parameter bound for this run
        with open(os.path.join(self.out_path, "param_bounds.json"), "w") as outfile: 
            json.dump(self.problem, outfile, indent=4)
    
    ############################################
    # Finalizing modules 
    ############################################

    def finalize(self):
        runtype = self.config['SALib']['method']
        
        if runtype == "Morris":
            self.plot_EET()
        else:
            print(f"Invalid runtype: {runtype}")
        
        # Add dotty plot module. Either here or in the above method
        # https://pynetlogo.readthedocs.io/en/latest/_docs/SALib_ipyparallel.html
                

    def plot_EET(self):
        """Plot elementary effects and std's for Morris analysis"""
        
        # Options for the graphic
        pltfont = {'fontname': 'DejaVu Sans', 'fontsize': 15}  # font for axes
        pltfont_leg = {'family': 'DejaVu Sans', 'size': 15}  # font for legend
        ms = 10  # Marker size
        col = np.array([[228, 26, 28], [55, 126, 184], [77, 175, 74],
                        [152, 78, 163], [255, 127, 0]]) / 256
        clrs = np.tile(col, (int(np.ceil(self.problem['num_vars'] / len(col))), 1))

        fig = plt.figure()

        # First plot EEs mean & std as circles:
        # Check the error bar definition
        for i in range(len(self.Si['mu_star'])):
            plt.errorbar(
                self.Si['mu_star'][i],     # x value
                self.Si['sigma'][i],                         # y value
                xerr=self.Si['mu_star_conf'][i],  # horizontal error (std deviation)
                fmt='o',                   # format for center marker
                markersize=ms, 
                color=clrs[i]
            )
            # plt.plot(
            #     self.Si['mu_star'][i], self.Si['sigma'][i], 'ok', markerfacecolor=clrs[i],
            #     markersize=ms, markeredgecolor='k'
            # )

        # Create legend:
        plt.legend(self.Si['names'], loc='best', prop=pltfont_leg)

        plt.xlabel('Mean of EEs', **pltfont)
        plt.ylabel('Standard deviation of EEs', **pltfont)
        plt.grid(linestyle='--')
        plt.xticks(**pltfont)
        plt.yticks(**pltfont)
        fig.set_size_inches(7, 7)
        plt.tight_layout()

        out_fn = 'EET.png'

        out_path_plot = os.path.join(self.out_path)
        plt.savefig(os.path.join(out_path_plot, out_fn), dpi=600, format='png')
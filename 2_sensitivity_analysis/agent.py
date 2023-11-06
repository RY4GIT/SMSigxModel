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
import datetime
import shutil


def read_SALib_config(config_SALib_path):
    """Reads configuration setting for SALib sensitivity analysis from text file."""

    # Read the data
    df = pd.read_csv(config_SALib_path)
    df_param_to_calibrate = df[df["calibrate"] == 1]

    # Convert it to dictionary
    config_SALib = {
        "num_vars": df_param_to_calibrate.shape[0],
        "names": df_param_to_calibrate["name"].tolist(),
        "bounds": df_param_to_calibrate[["lower_bound", "upper_bound"]].values.tolist(),
    }

    for i, bounds in enumerate(config_SALib["bounds"]):
        if bounds[1] <= bounds[0]:
            warnings.warn(
                f"The upper bound is smaller than the lower bound for the parameter {config_SALib['names'][i]}"
            )
    return config_SALib


class Agent_SALib_CFE:
    def __init__(self, config=None):
        """
        Initialize the sensitivity analysis using SALib for CFE model

        Sets up the initial state of the agent
        :param config:
        """

        self.config = config

        # Setting manual seed for reproducibility
        self.seed = 0

        # Define the parameter bounds for sensitivty analysis
        self.problem = read_SALib_config(config["PATHS"]["salib_config"])

        # Define output folder
        current_date = datetime.date.today().strftime("%Y-%m-%d")
        self.out_path = os.path.join(
            config["PATHS"]["homedir"],
            "results",
            f"{self.config['DATA']['site']}-{current_date}-{self.config['SALib']['like_measure']}-{self.config['SALib']['eval_variable']}",
        )
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)

    ############################################
    # Running sample & models
    ############################################

    def sample(self):
        """Run Morris sampling & analysis"""

        # Sample parameters
        N = int(self.config["Morris"]["N"])
        n_levels = int(self.config["Morris"]["n_levels"])
        self.sampled_params = morris_s.sample(
            self.problem, N=N, num_levels=n_levels, seed=self.seed
        )

        # Define number of runs
        self.nrun = N * (self.problem["num_vars"] + 1)
        print(
            f'Total runs: {self.nrun} \n Number of analyzed parameters: {self.problem["num_vars"]}\n'
        )

        return self.sampled_params

    def simulation(self, sampled_param_set):
        """
        Run sensitivity analysis, select the runner function based on the configuration

        Morris output:
            Si : A dictionary of sensitivity indices containing the following entries.
            mu : the mean elementary effect
            mu_star : the absolute of the mean elementary effect
            sigma : the standard deviation of the elementary effect
            mu_star_conf : the bootstrapped confidence interval
            names : the names of the parameters
        """

        self.model = CFEmodel(
            config=self.config, problem=self.problem, X=sampled_param_set
        )

        try:
            self.model.run()
            Y_i = self.model.evaluate()

            if np.isnan(Y_i):
                # Some parameters never impact soil moisture content: in this case replace Y_i with "perfect score"
                if self.config["SALib"]["eval_variable"] == "Soil Moisture Content":
                    if (self.config["SALib"]["like_measure"] == "KGE") or (
                        self.config["SALib"]["like_measure"] == "NashSutcliffe"
                    ):
                        return 1
                    elif "SeasonTrans" in self.config["SALib"]["like_measure"]:
                        return 0

            return Y_i

        except Exception as e:
            print(f"Error in simulation: {e}")
            return np.nan

    ############################################
    # Finalizing modules
    ############################################

    def finalize(self, results):
        """Finalizing the sensitivity analysis, select the finalizing function based on the configuration"""

        runtype = self.config["SALib"]["method"]

        if runtype == "Morris":
            self.analysis(results)
            self.plot_EET()
        else:
            print(f"Invalid runtype: {runtype}")

        self.remove_temp_files()
        # Add dotty plot module. Either here or in the above method
        # https://pynetlogo.readthedocs.io/en/latest/_docs/SALib_ipyparallel.html

    def remove_temp_files(self):
        directory = os.path.join(
            os.path.dirname(self.config["PATHS"]["cfe_config"]),
            "temporary_parameter_files_for_sensitivity_analysis",
        )
        shutil.rmtree(directory)

    def analysis(self, results):
        # Render the results returned from MP

        self.Y = np.array(results, dtype=float)

        # Analyze
        print("### Results ###")
        self.Si = morris_a.analyze(
            self.problem, self.sampled_params, self.Y, print_to_console=True
        )

        df_Si = self.Si.to_df()
        df_Si.to_csv(os.path.join(self.out_path, "Si.csv"))

        # Output the parameter bound for this run
        with open(os.path.join(self.out_path, "param_bounds.json"), "w") as outfile:
            json.dump(self.problem, outfile, indent=4)

    def plot_EET(self):
        """
        Plot elementary effects and std's for Morris analysis
            X-axis: mu_star : the absolute of the mean elementary effect
            Y-axis: sigma : the standard deviation of the elementary effect
        """

        # Options for the graphic
        pltfont = {"fontname": "DejaVu Sans", "fontsize": 15}  # font for axes
        pltfont_leg = {"family": "DejaVu Sans", "size": 15}  # font for legend
        ms = 10  # Marker size
        col = (
            np.array(
                [
                    [31, 120, 180],
                    [51, 160, 44],
                    [227, 26, 28],
                    [255, 127, 0],
                    [106, 61, 154],
                    [177, 89, 40],
                    [106, 61, 154],
                    [135, 135, 135],
                ]
            )
            / 256
        )
        clrs = np.tile(col, (int(np.ceil(self.problem["num_vars"] / len(col))), 1))

        fig = plt.figure()

        # First plot EEs mean & std as circles:
        # Check the error bar definition
        for i in range(len(self.Si["mu_star"])):
            if i <= 7:
                markertype = "ok"
            else:
                markertype = "x"
            plt.plot(
                self.Si["mu_star"][i],
                self.Si["sigma"][i],
                markertype,
                markerfacecolor=clrs[i],
                markersize=ms,
                markeredgecolor="k",
            )

        # Create legend:
        plt.legend(
            self.Si["names"],
            prop=pltfont_leg,
            loc="center left",
            bbox_to_anchor=(1, 0.5),
        )

        plt.xlabel("Mean of EEs", **pltfont)
        plt.ylabel("Standard deviation of EEs", **pltfont)
        plt.grid(linestyle="--")
        plt.xticks(**pltfont)
        plt.yticks(**pltfont)
        fig.set_size_inches(7, 7)
        plt.tight_layout()

        out_fn = "EET.png"

        out_path_plot = os.path.join(self.out_path)
        plt.savefig(os.path.join(out_path_plot, out_fn), dpi=600, format="png")

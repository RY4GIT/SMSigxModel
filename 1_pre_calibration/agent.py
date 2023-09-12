import os
import pandas as pd
import numpy as np
from model import CFEmodel
import datetime
import shutil
import spotpy


def read_spotpy_config(config_path):
    """Reads configuration setting for spotpy calibration analysis from text file."""

    # Read the data
    df = pd.read_csv(config_path)
    config_spotpy = df[df["calibrate"] == 1]

    return config_spotpy


def create_output_folder(home_dir, study_site):
    # Define output folder
    current_date = datetime.date.today().strftime("%Y-%m-%d")
    out_dir = os.path.join(
        home_dir,
        "results",
        f"{study_site}-{current_date}",
    )
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


class Spotpy_Agent:
    def __init__(self, config=None):
        """Initialize the Spotpy agent

        Args:
            config (_type_, optional): _description_. Defaults to None.
        """

        # Setup config files and output directories
        self.config = config
        self.nrun = int(self.config["DDS"]["N"])
        self.out_dir = create_output_folder(
            config["PATHS"]["homedir"], self.config["DATA"]["site"]
        )

        # Setting manual seed for reproducibility
        self.seed = 0
        np.random.seed(self.seed)

    def run(self):
        """Implement spotpy analysis"""

        # Setup spotpy
        spotpy_setup = Spotpy_setup(config=self.config)
        spotpy_runtype = self.config["spotpy"]["method"]
        if spotpy_runtype == "DDS":
            sampler = spotpy.algorithms.dds(
                spotpy_setup, dbname="raw_result_file", dbformat="ram"
            )
        else:
            print(f"Invalid runtype: {spotpy_runtype}")

        # https://github.com/thouska/spotpy/blob/c7f61c3333fd39e1f66b8211599279f80fe7fa6f/src/spotpy/describe.py#L37
        # https://github.com/thouska/spotpy/blob/c7f61c3333fd39e1f66b8211599279f80fe7fa6f/src/spotpy/algorithms/dds.py#L241
        sampler.sample(self.nrun)
        self.results = sampler.getdata()
        np.save(os.path.join(self.out_dir, "DDS_allresults.csv"), self.results)

    def finalize(self):
        self.reproduce_the_best_run(self.results)
        self.remove_temp_files()

    def reproduce_the_best_run(self, results):
        bestindex, bestobjf = spotpy.analyser.get_minlikeindex(results)
        best_model_run = results[bestindex]
        fields = [word for word in best_model_run.dtype.names if word.startswith("sim")]
        best_simulation = list(best_model_run[fields])
        np.save(os.path.join(self.out_dir, "DDS_bestrun.csv"), best_simulation)

    def remove_temp_files(self):
        directory = os.path.join(
            os.path.dirname(self.config["PATHS"]["cfe_config"]),
            "temporary_parameter_files_for_calibration",
        )
        shutil.rmtree(directory)


class Spotpy_setup:
    def __init__(self, config=None):
        self.config = config
        self.setup_params()

    def setup_params(self):
        param_bounds = read_spotpy_config(
            config_path=self.config["PATHS"]["spotpy_config"]
        )

        # setup calibration parameters
        self.params = [
            spotpy.parameter.Uniform(
                row["name"],
                low=row["lower_bound"],
                high=row["upper_bound"],
                optguess=row["optguess"],
            )
            for i, row in param_bounds.iterrows()
        ]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, x):
        self.model = CFEmodel(config=self.config, vector=x)

        try:
            self.model.run()
            self.sim_results = self.model.return_sim_runoff()
            return self.sim_results.values

        except Exception as e:
            print(f"Error in simulation: {e}")
            return np.nan  # TODO: fix the length?

    def evaluation(self, evaldates=False):
        vector_ = self.parameters()
        vector = dict()
        for item in vector_:
            vector[item[1]] = item[0]
        self.model = CFEmodel(config=self.config, vector=vector)
        self.obs_data = self.model.return_obs_runoff()
        return self.obs_data.values

    def objectivefunction(self, simulation, evaluation):
        if np.isnan(simulation.all()):
            self.obj_function = np.nan
        else:
            self.obj_function = spotpy.objectivefunctions.kge(
                evaluation[~np.isnan(evaluation)], simulation[~np.isnan(evaluation)]
            )
        return self.obj_function

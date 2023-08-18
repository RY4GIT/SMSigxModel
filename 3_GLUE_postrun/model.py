import os
import pandas as pd
pd.options.mode.chained_assignment = None
import json
import sys
sys.path.append(os.path.join(os.getcwd(), 'libs', 'py_cfe'))
from bmi_cfe import BMI_CFE

import shutil
from util import to_datetime


def duplicate_file(source_path, nth_run):
    """Duplicate a file appending nth_run to its name and return its path."""
    # Determine the directory and make a path for the new file
    directory = os.path.dirname(source_path)
    destination_dir = os.path.join(directory, "temporary_parameter_files_for_GLUE")
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)
    destination_path = os.path.join(destination_dir, f"config_cfe_{nth_run}.json")
    
    # Copy the source file to the new location
    shutil.copy2(source_path, destination_path)
    
    return destination_path  # Optionally return the path to the new file


class CFEmodel():
    def __init__(self, config=None):

        # Configs
        self.config = config


    def initialize(self, nth_run=0, sampled_params_set=None):
        """Initialize CFE model for the n-th GLUE run"""
        
        # Write the randomly-generated parameters to the config json file
        # Copy the CFE config file for sensitivty analysis 
        destination_path = duplicate_file(source_path = self.config['PATHS']['cfe_config'], nth_run=nth_run)

        # Get the model config file
        with open(destination_path) as data_file:
            self.cfe_cfg = json.load(data_file)

        # Overwrite the model config file
        for i in range(len(sampled_params_set)):
            if sampled_params_set[i][1] in ['bb', 'satdk', 'slop', 'satpsi', 'smcmax', 'wltsmc', 'D']:
                self.cfe_cfg["soil_params"][sampled_params_set[i][1]] = sampled_params_set[i][0]
            else:
                self.cfe_cfg[sampled_params_set[i][1]] = sampled_params_set[i][0]
                
        # Save the config file with new parameters
        with open(destination_path, 'w') as out_file:
            json.dump(self.cfe_cfg, out_file)
        
        # Initialize the CFE model with randomly sampled parameter
        self.cfe_instance = BMI_CFE(destination_path)
        self.cfe_instance.initialize()
            
    def run(self):
        """Run CFE model and return the synchronized timeseries of data & outputs"""
        self.cfe_instance.run_unit_test(plot=False, warm_up=True)
        
        # Returns the synchronized timeseries for the variables of interest
        obs = to_datetime(df=self.cfe_instance.unit_test_data, time_column="Time", format=self.config['DATA']['evaldata_time_format'])
        sim = to_datetime(df=self.cfe_instance.cfe_output_data, time_column="Time", format="%Y-%m-%d %H:%M:%S")
        df = pd.merge_asof(obs, sim, on = "Time")
        df = to_datetime(df, "Time", "%Y-%m-%d %H:%M:%S")

        obs_synced = pd.DataFrame()
        sim_synced = pd.DataFrame()
        
        self.var_names = ["Flow", "Soil Moisture Content", "Rainfall"]  
        for var_name in self.var_names:
            obs_synced[var_name] = df[var_name + "_x"].copy()
            sim_synced[var_name] = df[var_name + "_y"].copy()
        
        return obs_synced, sim_synced

    def get_observed_timeseries(self):
            # ==================================================
        # One more last run to get synced observed timeseries
        # ==================================================
        myCFE = BMI_CFE(self.config['PATHS']['cfe_config'])
        myCFE.initialize()
        
        obs = to_datetime(df=self.cfe_instance.unit_test_data, time_column="Time", format=self.config['DATA']['evaldata_time_format'])
        sim = to_datetime(df=self.cfe_instance.cfe_output_data, time_column="Time", format="%Y-%m-%d %H:%M:%S")
        obs_synced = pd.DataFrame()
        
        df = pd.merge_asof(obs, sim, on = "Time")
        df = to_datetime(df, "Time", "%Y-%m-%d %H:%M:%S")
        
        obs_synced = pd.DataFrame()
        sim_synced = pd.DataFrame()
        
        self.var_names = ["Flow", "Soil Moisture Content", "Rainfall"]  
        for var_name in self.var_names:
            obs_synced[var_name] = df[var_name + "_x"].copy()
            sim_synced[var_name] = df[var_name + "_y"].copy()
            
        self.df_obs_Q = pd.DataFrame(self.obs_synced["Flow"])
        self.df_obs_Q.set_axis(self.df_timeaxis, axis=0, inplace=True)
        self.df_obs_SM = pd.DataFrame(self.obs_synced["Soil Moisture Content"])
        self.df_obs_SM.set_axis(self.df_timeaxis, axis=0, inplace=True)
        
        return [self.df_obs_Q,  self.df_obs_SM]
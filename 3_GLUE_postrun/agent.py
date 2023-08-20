import os
import datetime
import ast

from model import CFEmodel
from GLUE import GLUE
from criteria import Criteria

# GLUE object
class Agent_GLUE_CFE_PostRun(object):
    def __init__(self, config=None):
        
        self.config=config
        criteria = Criteria(config=config)
        self.GLUE = GLUE(config=config, criteria=criteria)
        
        # Define the GLUE number of runs
        self.nrun = self.GLUE.nrun
        
        self.plot_figures = ast.literal_eval(self.config['GLUE']['plot_figures'])
            
    def evaluate(self):
        
        # Get the behavioral parameter sets using GLUE 
        behavioral_param_sets = self.GLUE.apply_criteria()
        
        # Plot figures 
        if self.plot_figures:
            self.GLUE.plot_parameter_distribution()
            self.GLUE.plot_parameter_dotty()
            self.GLUE.plot_parameter_interaction_dotty()
            
        return behavioral_param_sets

        
    def reproduce_a_behavioral_run(self, behavioral_params_set):
        
        # Preparations
        nth_run, behavioral_params_set = behavioral_params_set
        print(f'Processing {nth_run}/{self.nrun-1}')
        
        # Run CFE
        cfe_instance = CFEmodel(config=self.config)
        cfe_instance.initialize(nth_run=nth_run, behavioral_params_set=behavioral_params_set)
        _, sim = cfe_instance.run()

        # Return the variable of interest
        return [nth_run, sim['Flow'], sim['Soil Moisture Content']]


    def finalize(self, all_results):
        """Finalize"""
        
        # Check if there are behavioral runs 
        if len(all_results) == 0:
            print("There was no behavioral run")
            return 
        
        # Render the results from all GLUE runs to dataframe, store them in GLUE object
        self.GLUE.all_results_to_df(all_results)
        
        # Calculate the uncertainty bounds
        self.GLUE.get_uncertainty_bounds_for_all_variables()
        
        if self.plot_figures:
            
            # Get observed timeseries for a reference, use CFE BMI
            cfe_instance = CFEmodel(config=self.config)
            self.obs_Q, self.obs_SM = cfe_instance.get_observed_timeseries()
        
            # Plot the figure
            self.GLUE.plot_uncertainty_bounds(self.obs_Q, self.obs_SM)





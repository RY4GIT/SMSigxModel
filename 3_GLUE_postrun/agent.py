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
        
        # Define output folder
        # Get the current date in YYYY-MM-DD format
        current_date = datetime.date.today().strftime('%Y-%m-%d')
        self.out_path = os.path.join(config['PATHS']['homedir'], 'results', f"{self.config['DATA']['site']}-{current_date}", f"criteria_{criteria.id}")
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)
        
        # Define the GLUE number of runs
        self.nrun = self.GLUE.nrun
        
        self.plot_figures = ast.literal_eval(self.config['GLUE']['plot_figures'])
            
    def evaluate(self):
        
        # Get the behavioral parameter sets using GLUE 
        behavioral_param_sets = self.GLUE.apply_criteria()
        
        # Plot figures 
        if self.plot_figures:
            self.GLUE.plot_parameter_distribution()
            self.GLUE.plot_parametter_dotty()
            self.GLUE.plot_parameter_interaction_dotty()
            
        return behavioral_param_sets

        
    def reproduce_a_behavioral_run(self, behavioral_param):
        
        # Preparations
        nth_run, behavioral_params_set = behavioral_param
        print(f'Processing {nth_run}/{self.nrun-1}')
        
        # Run CFE
        cfe_instance = CFEmodel(config=self.config)
        cfe_instance.initialize(nth_run=nth_run, behavioral_params_set=behavioral_params_set)
        obs, sim = cfe_instance.run()

        # Get the variable of interest
        self.sim_behavioral_Q = sim['Flow']
        self.sim_behavioral_SM = sim['Soil Moisture Content']
        
        return [nth_run, self.df_behavioral_Q, self.df_behavioral_SM]


    def finalize(self, all_results):
        
        # Render the results from all GLUE runs to dataframe, store them in GLUE object
        self.GLUE.all_results_to_df(all_results)
        
        # Calculate the uncertainty bounds
        self.GLUE.calc_uncertainty_bounds()
        
        if self.plot_figures:
            
            # Get observed timeseries for a reference, use CFE BMI
            cfe_instance = CFEmodel(config=self.config)
            self.df_obs_Q, self.df_obs_SM = cfe_instance.get_observed_timeseries()
        
            # Plot the figure
            self.GLUE.plot_uncertainty_bounds(self.df_obs_Q, self.df_obs_SM)





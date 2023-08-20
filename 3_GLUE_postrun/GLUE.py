
import pandas as pd
import numpy as np
import os
from util import triangle_weight, weighted_quantile
import matplotlib.pyplot as plt
import datetime

class GLUE(object):
    def __init__(self, config=None, criteria=None):
        
        self.config = config
        
        ################################
        # Get configuration for GLUE analysis
        ################################
        
        # Quantile
        self.quantiles = [float(x.strip()) for x in self.config.get('GLUE', 'quantiles').split(',')]
        
        # Get GLUE criteria
        self.criteria = criteria
    
        # Define output folder
        # Get the current date in YYYY-MM-DD format
        current_date = datetime.date.today().strftime('%Y-%m-%d')
        self.out_path = os.path.join(config['PATHS']['homedir'], 'results', f"{self.config['DATA']['site']}-{current_date}", f"criteria_{criteria.id}")
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)
            
        # Keyword arguments for to_csv
        self.kwargs_dict_csv = {
            "sep": ',',
            "header": True,
            "index": True,
            "encoding": 'utf-8',
            "na_rep": 'nan'
        }
            
        ################################
        # Get configuration used in GLUE prerun
        ################################
        
        # Get the parameter sets and evaluation metrics from GLUE pre-run
        GLUEprerun_output_path = self.config['PATHS']['GLUE_prerun_output_path']
        self.prior_params = pd.read_csv(os.path.join(GLUEprerun_output_path, 'prior_parameters.csv'), index_col=0)
        self.prior_eval_metrics = pd.read_csv(os.path.join(GLUEprerun_output_path, 'evaluation_metrics.csv'), index_col=0)
        self.prior_eval_metrics_monthly = pd.read_csv(os.path.join(GLUEprerun_output_path, 'evaluation_metrics_monthly.csv'), index_col=0)
        
        # Number of runs
        self.nrun = len(self.prior_params)

    
    def judge_behavioral(self, metric_values, threshold, behavioral_logical_operation):
        if behavioral_logical_operation == "metric value is more than threshold":
            return metric_values > threshold
        elif behavioral_logical_operation == "metric value is less than threshold":
            return metric_values < threshold
    
    def apply_criteria(self):
        """Evaluate GLUE prerun reuslts based on the criteria"""

        #########################
        # Main GLUE procedure 
        #########################
        
        # Judget behavioral vs. non-behavioral for all criterion
        for _, criterion in self.criteria.full_criteria.items():
            self.prior_eval_metrics[criterion['metrics_fullname'] + '_Behavioral'] = self.judge_behavioral(
                metric_values=self.prior_eval_metrics[criterion['metrics_fullname']], 
                threshold=criterion['threshold'], 
                behavioral_logical_operation=criterion['operation']
            )
            
        # Judget behavioral vs. non-behavioral based on all criteria (multi-criteria)
        glue_judge_column = self.prior_eval_metrics.filter(like='_Behavioral').columns
        self.prior_eval_metrics['GLUE_results_Behavioral'] = self.prior_eval_metrics[glue_judge_column].all(axis='columns')
        
        # Get the run_id of behavioral runs
        self.behavioral_run_ids = self.prior_eval_metrics[self.prior_eval_metrics['GLUE_results_Behavioral']].index
        
        #########################
        # Get the behavioral run attributes and save
        #########################
        
        # Get evaluation metrics from only behavioral runs 
        self.posterior_eval_metrics = self.prior_eval_metrics.iloc[self.behavioral_run_ids].copy()
        self.posterior_eval_metrics_mo = self.prior_eval_metrics_monthly[self.prior_eval_metrics_monthly.index.isin(self.behavioral_run_ids.values)].copy()
        
        # Save 
        self.posterior_eval_metrics.to_csv(os.path.join(self.out_path, 'post_evaluations.csv'), **self.kwargs_dict_csv)
        self.posterior_eval_metrics_mo.to_csv(os.path.join(self.out_path, 'post_evaluations_monthly_metrics.csv'), **self.kwargs_dict_csv)
        
        # Get only behavioral parameters 
        self.posterior_params = self.prior_params.iloc[self.behavioral_run_ids].copy()

        # Get the bahvioral parmaeter sets as a list
        behavioral_params = [[run_id, self.posterior_params.loc[run_id]] for run_id in self.behavioral_run_ids]
        
        return behavioral_params
        
    def all_results_to_df(self, all_results):
        """Render all reuslts to dataframe"""
        
        self.run_id_behavioral = [result[0] for result in all_results]

        for i, result in enumerate(all_results):
            if not hasattr(self, 'behavioral_Q'):
                self.behavioral_Q = result[1]
                self.behavioral_SM = result[2]
            else:
                self.behavioral_Q = pd.concat([self.behavioral_Q, result[1]], axis=1)
                self.behavioral_SM = pd.concat([self.behavioral_SM, result[2]], axis=1)

        self.behavioral_Q.set_axis(self.run_id_behavioral, axis=1, inplace=True)
        self.behavioral_SM.set_axis(self.run_id_behavioral, axis=1, inplace=True)

    def calc_weights(self, eval_metrics_values, threshold, metrics_fullname):
        if 'SeasonTrans' in metrics_fullname:
            weight = triangle_weight(eval_metrics_values,
                                                a=-1 * threshold, b=0,
                                                c=threshold)
        else:
            weight = ((eval_metrics_values - threshold) / 
                            sum(eval_metrics_values - threshold)).to_numpy()
        return weight
        
    def get_weights_for_multi_criteria(self):
        # Weight matrix: [number of behavioral runs] - by - [number of criterial]
        weights = np.empty((len(self.behavioral_run_ids), len(self.criteria.full_criteria)))

        # Get weights for each criteria
        for i, (_, criterion) in enumerate(self.criteria.full_criteria.items()):
            weights[:, i] = self.calc_weights(self.posterior_eval_metrics[criterion['metrics_fullname']], criterion['threshold'], criterion['metrics_fullname'])
            
        # Average weight for all criteria
        return np.mean(weights, axis=1)
    
    def calc_uncertainty_bounds(self, var_attr, avg_weight):
        
        # Get the dataseries based on the variable attribute names
        df_behavioral = getattr(self, f'behavioral_{var_attr}').copy()
        np_behavioral = df_behavioral.to_numpy(copy=True)

        # Get the quantiles of dataseries for each timstep
        _dataseries_quantiles = np.array([
            weighted_quantile(values=row, quantiles=self.quantiles, sample_weight=avg_weight) 
            for row in np_behavioral
        ])
        
        # Render them as pandas dataframe
        dataseries_quantiles = pd.DataFrame(_dataseries_quantiles, index=df_behavioral.index, columns=['lowerlim', 'median', 'upperlim'])
        setattr(self, f"quantile_{var_attr}", dataseries_quantiles)

        # Save results
        if hasattr(self, f"quantile_{var_attr}"):
            dataseries_quantiles.to_csv(os.path.join(self.out_path, f'quantiles_{var_attr}.csv'), **self.kwargs_dict_csv)
        
    def get_uncertainty_bounds_for_all_variables(self):
        
        # Calculate average weight for multi-criteria
        avg_weight = self.get_weights_for_multi_criteria()
    
        variable_map = {
            "Flow": "Q",
            "Soil Moisture Content": "SM"
        }

        # Calculate uncertainty bounds for each variables of interest
        for _, var_attr in variable_map.items():
            self.calc_uncertainty_bounds(var_attr, avg_weight)
    
    def plot_parameter_distribution(self):
        
        nparas = len(self.prior_params.columns)
        
        # Setup figure
        f, axes = plt.subplots(nparas, 1, figsize=(4 * 4, 4 * nparas))
        
        for i, target_para in enumerate(self.prior_params.columns):
            ax = axes[i]
            
            # Plot histograms
            self.prior_params[target_para].hist(bins=10, alpha=0.4, ax=ax, color="#3182bd", label="Prior")
            
            if hasattr(self, 'posterior_params'):
                self.posterior_params[target_para].hist(bins=10, alpha=0.8, ax=ax, color="#3182bd", label="Posterior")
            
            ax.set_xlabel(target_para)
            ax.legend()
            
            if i != 0:
                ax.yaxis.set_visible(False)

        # Save figure
        f.savefig(os.path.join(self.out_path, 'param_dist.png'), dpi=600)
        plt.close(f)

    def plot_parameter_dotty(self):

        nparas = len(self.prior_params.columns)
        eval_metrics = [col for col in self.posterior_eval_metrics.columns if '_Behavioral' not in col]

        for target_eval in eval_metrics:
            f, axes = plt.subplots(nparas, 1, figsize=(4 * 4, 4 * nparas), constrained_layout=True)

            for i, target_para in enumerate(self.prior_params.columns):
                ax = axes[i]
                ax.scatter(self.posterior_params[target_para], self.posterior_eval_metrics[target_eval], alpha=0.5)
                ax.tick_params(axis='both', which='major', labelsize=16)
                ax.set_xlabel(target_para, fontsize=16)
                ax.set_ylabel(target_eval, fontsize=16)

            save_path = os.path.join(self.out_path, f"param_dotty_{target_eval}.png")
            f.savefig(save_path, dpi=600)
            plt.close(f)
                
    def plot_parameter_interaction_dotty(self):
        # List of parameters of interest
        param_interest = self.posterior_params.columns.to_list()

        n_params = len(param_interest)
        
        # Create a grid of subplots
        f, axes = plt.subplots(n_params, n_params, figsize=(20, 20), constrained_layout=True)

        for i, para1 in enumerate(param_interest):
            for j, para0 in enumerate(param_interest):
                ax = axes[i, j]
                ax.scatter(self.posterior_params[para0], self.posterior_params[para1], alpha=0.5)
                
                if i == 0:
                    ax.xaxis.set_label_position('top')
                    ax.set_xlabel(para0)
                if j == 0:
                    ax.set_ylabel(para1)
                ax.tick_params(direction="in")

        save_path = os.path.join(self.out_path, "param_dotty_interaction.png")
        f.savefig(save_path, dpi=600)
        plt.close(f)

        
    def plot_uncertainty_bounds(self, observed_Q, observed_SM):

        settings = {
            "Flow": {
                "df_simrange": self.quantile_Q,
                "df_obs": observed_Q,
                "obs_label": 'Observed flow',
                "y_label": 'Total flow [mm/hour]',
                "title": 'Total flow for behavioral runs',
                "fn": 'quantiles_Q.png',
                "yscale": 'log',
                "ylim": [1E-07, 1E-01]
            },
            "Soil Moisture Content": {
                "df_simrange": self.quantile_SM,
                "df_obs": observed_SM,
                "obs_label": 'Observed soil moisture',
                "y_label": 'Volumetric Soil Moisture Content [m^3/m^3]',
                "title": 'Soil moisture for behavioral runs',
                "fn": 'quantiles_SM.png'
            }
        }

        for _, setting in settings.items():
            f2, ax2 = plt.subplots(figsize=(8, 6))
            
            setting["df_simrange"]['lowerlim'].plot(color='black', alpha=0.2, ax=ax2, label=f'{self.quantiles[0]*100} percentile')
            setting["df_simrange"]['upperlim'].plot(color='black', alpha=0.2, ax=ax2, label=f'{self.quantiles[2]*100} percentile')
            setting["df_obs"].plot(color='black', alpha=1, ax=ax2, label=setting["obs_label"])

            plt.fill_between(setting["df_simrange"].index, setting["df_simrange"]['upperlim'], setting["df_simrange"]['lowerlim'],
                            facecolor='green', alpha=0.2, interpolate=True, label='Predicted range')

            if "yscale" in setting:
                ax2.set_yscale(setting["yscale"])
                ax2.set_ylim(setting["ylim"])

            ax2.set_xlabel('Time')
            ax2.set_ylabel(setting["y_label"])
            ax2.set_title(setting["title"])
            ax2.legend()

            f2.autofmt_xdate()
            f2.savefig(os.path.join(self.out_path, setting["fn"]), dpi=600)
            plt.close(f2)
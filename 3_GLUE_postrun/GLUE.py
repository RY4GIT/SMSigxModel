
import pandas as pd
import numpy as np
import os
from util import triangle_weight, weighted_quantile
import matplotlib.pyplot as plt
import datetime

class GLUE(object):
    def __init__(self, config=None, criteria=None):
        
        self.config = config
        
        # Get GLUE criteria
        self.criteria = criteria
        
        # Get previous GLUE configuration
        
        # Get the parameter sets and evaluation metrics from GLUE pre-run
        GLUEprerun_output_path = self.config['PATHS']['GLUE_prerun_output_path']
        self.prior_params = pd.read_csv(os.path.join(GLUEprerun_output_path, 'prior_parameters.csv'), index_col=0)
        self.prior_eval_metrics = pd.read_csv(os.path.join(GLUEprerun_output_path, 'evaluation_metrics.csv'), index_col=0)
        self.prior_eval_metrics_monthly = pd.read_csv(os.path.join(GLUEprerun_output_path, 'evaluation_metrics_monthly.csv'), index_col=0)
        
        # Number of runs
        self.nrun = len(self.prior_params)
        
        # Quantile
        self.quantiles = [float(x.strip()) for x in self.config.get('GLUE', 'quantiles').split(',')]
        
        # Define output folder
        # Get the current date in YYYY-MM-DD format
        current_date = datetime.date.today().strftime('%Y-%m-%d')
        self.out_path = os.path.join(config['PATHS']['homedir'], 'results', f"{self.config['DATA']['site']}-{current_date}", f"criteria_{criteria.id}")
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)
    
    def judge_behavioral(self, dataseries, threshold, behavioral_logical_operation):
        if behavioral_logical_operation == "metric value is more than threshold":
            return dataseries > threshold
        elif behavioral_logical_operation == "metric value is less than threshold":
            return dataseries < threshold
    
    def apply_criteria(self):
        """Evaluate GLUE prerun reuslts based on the criteria"""

        #########################
        # Main GLUE procedure 
        #########################
        
        # Judget behavioral vs. non-behavioral by each criterion
        for _, criterion in self.criteria.full_criteria.items():
            self.prior_eval_metrics[criterion['metrics_fullname'] + '_Behavioral'] = self.judge_behavioral(
                dataseries=self.prior_eval_metrics[criterion['metrics_fullname']], 
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
        self.posterior_eval_metrics.to_csv(os.path.join(self.out_path, 'post_evaluations.csv'), sep=',', header=True,
                                    index=True, encoding='utf-8', na_rep='nan')
        self.posterior_eval_metrics_mo.to_csv(os.path.join(self.out_path, 'post_evaluations_monthly_metrics.csv'), sep=',', header=True,
                                    index=True, encoding='utf-8', na_rep='nan')
        
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

    def calc_weights(self):
        weights = np.empty((len(self.posterior_eval_metrics), len(self.eval_names)))

        j = 0
        for i, criteria in enumerate(self.eval_criteria):
            metric = criteria['metric']
            if metric in ["NSE", "KGE"]:
                weights[:, j] = ((self.posterior_eval_metrics[self.eval_names[j]] - criteria['threshold']) / 
                                sum(self.posterior_eval_metrics[self.eval_names[j]] - criteria['threshold'])).to_numpy()
                j += 1
            elif metric == "season_transition":
                for k in range(4):
                    weights[:, j + k] = triangle_weight(self.posterior_eval_metrics[self.eval_names[j + k]],
                                                        a=-1 * criteria['threshold'], b=0,
                                                        c=criteria['threshold'])
                j += 4
        return np.mean(weights, axis=1)
    
    def calc_uncertainty_bounds(self, plot=True):
        
        avg_weight = self.calc_weights()
    
        variable_map = {
            "Flow": "df_Q_simrange",
            "Soil Moisture Content": "df_SM_simrange"
        }

        for var_name, var_attr in variable_map.items():
            df_behavioral = getattr(self, f'df_behavioral_{var_name[0]}').copy()
            np_behavioral = df_behavioral.to_numpy(copy=True)

            quantile = np.array([
                weighted_quantile(values=row, quantiles=self.quantiles, sample_weight=avg_weight) 
                for row in np_behavioral
            ])
            
            df_simrange = pd.DataFrame(quantile, index=df_behavioral.index, columns=['lowerlim', 'median', 'upperlim'])
            setattr(self, var_attr, df_simrange)

            if hasattr(self, var_attr):
                df_simrange.to_csv(os.path.join(self.out_path, f'quantiles_{var_name[0]}.csv'), sep=',', 
                                header=True, index=True, encoding='utf-8', na_rep='nan')
    
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
                "df_simrange": self.df_Q_simrange,
                "df_obs": self.df_obs_Q,
                "obs_label": 'Observed flow',
                "y_label": 'Total flow [mm/hour]',
                "title": 'Total flow for behavioral runs',
                "fn": 'Q_range.png',
                "yscale": 'log',
                "ylim": [1E-07, 1E-01]
            },
            "Soil Moisture Content": {
                "df_simrange": self.df_SM_simrange,
                "df_obs": self.df_obs_SM,
                "obs_label": 'Observed soil moisture',
                "y_label": 'Volumetric Soil Moisture Content [m^3/m^3]',
                "title": 'Soil moisture for behavioral runs',
                "fn": 'SM_range.png'
            }
        }

        for var_name, setting in settings.items():
            f2, ax2 = plt.subplots(figsize=(8, 6))
            
            setting["df_simrange"]['lowerlim'].plot(color='black', alpha=0.2, ax=ax2, label=f'{quantiles[0]*100} percentile')
            setting["df_simrange"]['upperlim'].plot(color='black', alpha=0.2, ax=ax2, label=f'{quantiles[2]*100} percentile')
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
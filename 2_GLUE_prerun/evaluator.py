import os
import sys

sys.path.append(os.path.join(os.getcwd(), 'libs', 'SMSig'))
from sig_seasontrans import SMSig

import numpy as np
import pandas as pd
from statistics import median  

import spotpy    

class Evaluator():
    
    def __init__(self, simulation=None, observation=None) -> None:
        self.flow_var_name = "Flow"
        self.soilmoisture_var_name = "Soil Moisture Content"
        
        self.observation = observation
        self.simulation = simulation
        
        
    def evaluate(self):
        """Evaluate simulation results."""
        
        eval_hourly = self._evaluate_hourly()
        eval_monthly = self._evaluate_monthly()
        
        return eval_hourly, eval_monthly

    def _evaluate_hourly(self):
        """Evaluate metrics at model (hourly) timestep."""
        
        metrics = {
            "NSE on Flow": self.calc_NSE(self.flow_var_name),
            "NSE on Soil": self.calc_NSE(self.soilmoisture_var_name),
            "KGE on Flow": self.calc_KGE(self.flow_var_name),
            "KGE on Soil": self.calc_KGE(self.soilmoisture_var_name)
        }
        
        season_trans = self.calc_SeasonTrans(self.soilmoisture_var_name)
        metrics.update({
            "SeasonTrans of Soil dry2wet_start": season_trans[0],
            "SeasonTrans of Soil dry2wet_end": season_trans[1],
            "SeasonTrans of Soil wet2dry_start": season_trans[2],
            "SeasonTrans of Soil wet2dry_end": season_trans[3]
        })

        return pd.DataFrame(metrics, index=[0])

    def _evaluate_monthly(self):
        """Evaluate metrics at model (monthly) timestep."""
        
        self.observation_monthly = self.df_to_monthly_timestep(self.observation)
        self.simulation_monthly = self.df_to_monthly_timestep(self.simulation)
        
        median_flow_obs = median(self.observation['Flow'])
        high_flow_obs = 9 * median_flow_obs
        
        metrics = {
            'Q_mean_obs': self.calc_monthly_Q_mean(self.flow_var_name, self.observation_monthly),
            'Q_mean_sim': self.calc_monthly_Q_mean(self.flow_var_name, self.simulation_monthly),
            'high_flow_freq_obs': self.calc_monthly_high_flow_freq(self.flow_var_name, self.observation, high_flow_obs),
            'high_flow_freq_sim': self.calc_monthly_high_flow_freq(self.flow_var_name, self.simulation, high_flow_obs),
            'RR_obs': self.calc_monthly_RR(self.flow_var_name, self.observation_monthly, self.observation['Rainfall']),
            'RR_sim': self.calc_monthly_RR(self.flow_var_name, self.simulation_monthly, self.observation['Rainfall'])
        }
        
        return pd.concat(metrics.values(), axis=1, keys=metrics.keys())
        
    def calc_NSE(self, variable):
        return spotpy.objectivefunctions.nashsutcliffe(self.observation[variable], self.simulation[variable])
        
    def calc_KGE(self, variable):
        return spotpy.objectivefunctions.kge(self.observation[variable], self.simulation[variable])
        
    def calc_SeasonTrans(self, variable):
        
        # Calculate metrics for observed timeseries
        sig_obs = SMSig(
                ts_time= self.observation.index.to_numpy(),
                ts_value= self.observation[variable].to_numpy(),
                plot_results=False,
                plot_label="obs"
            )
        sig_obs.movmean()
        t_valley = sig_obs.calc_sinecurve()
        season_trans_obs, _, _ = sig_obs.calc_seasontrans(t_valley=t_valley)

        # Calculate metrics for SIMULATED timeseries
        sig_sim = SMSig(
                ts_time=self.simulation.index.to_numpy(),
                ts_value=self.simulation[variable].to_numpy(),
                plot_results=False,
                plot_label="sim"
            )
        sig_sim.movmean()
        season_trans_sim, _, _ = sig_sim.calc_seasontrans(t_valley=t_valley)

        # Get the deviations in seasonal transition dates between simulated and observed timeseries
        diff = season_trans_sim - season_trans_obs
        abs_diff_SeasonTransDate = abs(np.nanmean(diff, axis=0))

        return abs_diff_SeasonTransDate
        
    def df_to_monthly_timestep(self, df):
        df_monthly = df.resample('M').sum().copy()
        df_monthly.drop(df_monthly.tail(1).index, inplace=True)
        df_monthly.drop(df_monthly.head(1).index, inplace=True)
        return df_monthly
        
    def calc_monthly_Q_mean(self, variable, df):
        return df[variable].resample('M').mean()

    def calc_monthly_high_flow_freq(self, variable, df, high_flow_criteria):
        monthly_high_flow_freq = df[variable].resample('M').apply(lambda x: sum(x > high_flow_criteria) / len(x))
        return monthly_high_flow_freq

    def calc_monthly_RR(self, variable, df, precip):
        monthly_flow = df[variable].resample('M').sum()
        monthly_precip = precip.resample('M').sum()
        return monthly_flow / monthly_precip
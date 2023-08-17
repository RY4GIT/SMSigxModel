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
        self.soilmoisture_var_name = "Soil Moisture Content"
        self.flow_var_name = "Flow"
                
        self.simulation = simulation
        self.observation = observation
        
    def evaluate(self):
        """Evaluate simulation results"""
        
        ################################################
        # Evaluation metrics at model (hourly) timestep
        ################################################
        
        NSE_Flow = self.calc_NSE(variable=self.flow_var_name)
        NSE_Soil = self.calc_NSE(variable=self.soilmoisture_var_name)
        KGE_Flow = self.calc_KGE(variable=self.flow_var_name)
        KGE_Soil = self.calc_KGE(variable=self.soilmoisture_var_name)
        SeasonTrans_Soil = self.calc_SeasonTrans(variable=self.soilmoisture_var_name)
        
        eval_hourly = pd.DataFrame({
            "NSE on Flow": NSE_Flow,
            "NSE on Soil": NSE_Soil,
            "KGE on Flow": KGE_Flow,
            "KGE on Soil": KGE_Soil,
            "SeasonTrans_Soil dry2wet_start": SeasonTrans_Soil[0],
            "SeasonTrans_Soil dry2wet_end": SeasonTrans_Soil[1],
            "SeasonTrans_Soil wet2dry_start": SeasonTrans_Soil[2],
            "SeasonTrans_Soil wet2dry_end": SeasonTrans_Soil[3]
        })
        
        ################################################
        # Evaluation metrics at model (hourly) timestep
        ################################################
        
        self.simulation_monthly = self.df_to_monthly_timestep(self.simulation)
        self.observation_monthly = self.df_to_monthly_timestep(self.observation)
        
        # Calculate mean soil moisture
        Q_mean_obs = self.calc_monthly_Q_mean(variable=self.flow_var_name, df=self.observation_monthly)
        Q_mean_sim = self.calc_monthly_Q_mean(variable=self.flow_var_name, df=self.simulation_monthly)
        
        median_flow_obs = median(self.observation['Flow'])
        high_flow_obs = 9 * median_flow_obs
        
        high_flow_freq_obs = self.calc_monthly_high_flow_freq(variable=self.flow_var_name, df=self.observation_monthly, high_flow_criteria=high_flow_obs)
        high_flow_freq_sim = self.calc_monthly_high_flow_freq(variable=self.flow_var_name, df=self.simulation_monthly, high_flow_criteria=high_flow_obs)
        
        RR_obs = self.calc_monthly_RR(variable=self.flow_var_name, df=self.observation_monthly, precip=self.observation['Rainfall'])
        RR_sim = self.calc_monthly_RR(variable=self.flow_var_name, df=self.simulation_monthly, precip=self.observation['Rainfall'])
        
        eval_monthly = [Q_mean_obs, Q_mean_sim, high_flow_freq_obs, high_flow_freq_sim, RR_obs, RR_sim]
        
        return eval_hourly, eval_monthly
        
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
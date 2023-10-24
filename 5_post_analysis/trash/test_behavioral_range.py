# A main module to run various analysis with CFE model
import os
import sys

os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/5_GLUE_model")

# Import libraries
import spotpy
import numpy as np
import pandas as pd
import shutil
import json

sys.path.append("../libs/cfe/cfe_py")
import cfe
import bmi_cfe

sys.path.append("../libs/SMSig")
from sig_seasontrans import SMSig

sys.path.append("../libs/")
from glue_cfe_mp import MyGLUE

out_path = (
    r"G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\6_out\Mahurangi\ex111"
)

# Read the behavioral runs
df_SM = pd.read_csv(os.path.join(out_path, "simulated_SM.csv"), index_col=0)
df_Q = pd.read_csv(os.path.join(out_path, "simulated_Q.csv"), index_col=0)
df_glue_results = pd.read_csv(os.path.join(out_path, "glue_results.csv"), index_col=0)
df_eval = pd.read_csv(os.path.join(out_path, "evaluations.csv"), index_col=0)
behavioral_run_id_index = df_glue_results.index[df_glue_results["Behavioral"].values]
df_behavioral_SM = df_SM[behavioral_run_id_index.astype(str)].copy()
df_behavioral_Q = df_Q[behavioral_run_id_index.astype(str)].copy()
df_post_eval = df_eval.iloc[behavioral_run_id_index]

# Create GLUE instance
# Create a GLUE instance
out_path_2 = os.path.join(out_path, "quantile_25_75")
if not os.path.exists(out_path_2):
    os.mkdir(out_path_2)
glue_instance = MyGLUE(out_path=out_path_2)

# Overwrite the behavioral runs
glue_instance.df_post_eval = df_post_eval
glue_instance.eval_names = ["NSE on Flow"]
glue_instance.eval_criteria = {
    0: {"variable_to_analyze": "Flow", "metric": "NSE", "threshold": 0.5}
}
glue_instance.df_behavioral_Q = df_behavioral_Q
glue_instance.df_behavioral_SM = df_behavioral_SM
glue_instance.df_obs_Q = df_Q["0"]
glue_instance.df_obs_Q.rename({"0": "Flow"})
glue_instance.df_obs_SM = df_SM["0"]

# Try postporcess
glue_instance.post_process()

# Output the plot
glue_instance.plot("timeseries")

# A main module to run various analysis with CFE model

# Import libraries
import os
import sys
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/5_GLUE_model")
pool_nprocess = 3
os.environ['NUMEXPR_MAX_THREADS'] = '3'
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import cProfile, pstats, io # https://docs.python.org/3/library/profile.html#module-cProfile
from pstats import SortKey
import multiprocessing as mp
import snakeviz
import spotpy
import pandas as pd
sys.path.append("../libs/")
from glue_post_analysis import MyGLUEPost

def main(out_path='', config_path_CFE='', path_GLUE_output='', eval_criteria=dict()):
    # To reproduce behavioral runs and evaluate
    print("Post evaluation of GLUE analysis")

    # Create output directly if not exist
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # Create a GLUE instance
    glue_instance = MyGLUEPost(
        out_path=out_path,
        config_path_CFE=config_path_CFE,
        path_GLUE_output=path_GLUE_output,
        eval_criteria=eval_criteria
    )
    
    behavioral_params = glue_instance.evaluation(plot=False)
    

    # Start multiple runs in multiprocessing
    # pool = mp.Pool(processes=pool_nprocess)
    # all_results = pool.map(glue_instance.reproduce_behavioral_run, behavioral_params)
    # # all_results = pool.map(glue_instance.simulation, behavioral_params[0:3])
    # pool.close()
    # pool.join()

        # glue_instance.post_process()
    # Post-process the results
    # glue_instance.calc_uncertainty_bounds(all_results)
    
    # glue_instance.post_eval(eval_criteria_1)
    

    
    print(f'Finished the run')
    print(f'Saved results to: {out_path}')

if __name__ == '__main__':

    eval_criteria_1 = {
        0: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'KGE', 'threshold': 0.5},
        1: {'variable_to_analyze': 'Flow', 'metric': 'KGE', 'threshold':  0.5},
        2: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'season_transition', 'threshold': 30}
    }
        
    main(
        out_path= r'../8_out/Mahurangi/ex111',
        config_path_CFE= r'../2_data_input/Mahurangi/parameters/config_cfe_template.json',
        path_GLUE_output= r'../6_out/Mahurangi/ex111',
        eval_criteria=eval_criteria_1
    )

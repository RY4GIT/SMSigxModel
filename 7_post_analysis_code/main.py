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

def main(out_path='', config_path_CFE='', path_GLUE_output=''):
    # To reproduce behavioral runs and evaluate

    # Create output directly if not exist
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    eval_criterion_1 = {
        0: {'variable_to_analyze': 'Flow', 'metric': 'KGE', 'threshold': 0.5}
        }
    
    eval_criterion_2 = {
        0: {'variable_to_analyze': 'Flow', 'metric': 'KGE', 'threshold': 0.5},
        1: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'KGE', 'threshold':  0.5}
        }
    
    eval_criterion_3 = {
        0: {'variable_to_analyze': 'Flow', 'metric': 'KGE', 'threshold': 0.5},
        1: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'season_transition', 'threshold': 30}
    }
    
    eval_criteria = [eval_criterion_1, eval_criterion_2, eval_criterion_3]
    senario_ids = [1, 2, 3]
        
    # Create a GLUE instance
    for i in range(len(senario_ids)):
        
        eval_criterion = eval_criteria[i]
        senario_id = senario_ids[i]
        print(f"Post evaluation of GLUE analysis: Criterion {senario_id}")
            
        glue_instance = MyGLUEPost(
            out_path=out_path,
            config_path_CFE=config_path_CFE,
            path_GLUE_output=path_GLUE_output,
            eval_criteria=eval_criterion,
            senario_id=senario_id
        )
        
        behavioral_params = glue_instance.post_evaluation(plot=True)
        
        # Start multiple runs in multiprocessing
        print('--- Reproducing behavioral runs ---')
        pool = mp.Pool(processes=pool_nprocess)
        all_results = pool.map(glue_instance.reproduce_behavioral_run, behavioral_params)
        # all_results = pool.map(glue_instance.simulation, behavioral_params[0:3])
        pool.close()
        pool.join()

        # Calculate the uncertainty bounds
        glue_instance.calc_uncertainty_bounds(all_results, plot=True)
    
    print(f'Finished the run')
    print(f'Saved results to: {glue_instance.out_path_per_senario}')

if __name__ == '__main__':
        
    main(
        out_path= r'../8_out/Mahurangi/ex111',
        config_path_CFE= r'../2_data_input/Mahurangi/parameters/config_cfe_template.json',
        path_GLUE_output= r'../6_out/Mahurangi/ex111'
    )

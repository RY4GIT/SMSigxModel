# A main module to run various analysis with CFE model

# Import libraries
import os
import configparser
from agent import Agent_GLUE_CFE_PostRun
import multiprocessing as mp

def read_config(config_path): 
    config = configparser.ConfigParser()
    config.read(config_path)
    return config

def main(config_path):

    print(f'#### Start GLUE experiment - Postrun ###')
    
    # Read config
    config = read_config(config_path=config_path)
    
    # Create a GLUE-CFE agent
    glue_experiment = Agent_GLUE_CFE_PostRun(config)
    behavioral_param_sets = glue_experiment.evaluate()
    
    # Start multiple runs in multiprocessing
    print('--- Reproducing behavioral runs ---')
    
    # Just to excute a single model run for a test
    if glue_experiment.nrun > 1: 
        print('--- Multithreadding mode ---')
        
        # MP settings
        print("Number of cpu : ", mp.cpu_count())
        pool_nprocess = int(config['Multiprocessing']['pool_nprocess'])
        os.environ['NUMEXPR_MAX_THREADS'] = config['Multiprocessing']['NUMEXPR_MAX_THREADS']
        
        # Run MP
        print('--- Reproducing behavioral runs started ---')
        pool = mp.Pool(processes=pool_nprocess)
        all_results = pool.map(glue_experiment.reproduce_a_behavioral_run, behavioral_param_sets)
        pool.close()
        pool.join()
        print(f"--- Finished reproducing runs ---")
        
    elif glue_experiment.nrun == 1:
        print('--- Single run mode ---')
        all_results = [glue_experiment.reproduce_a_behavioral_run(behavioral_param_set=behavioral_param_sets[0])]
    
    # To execute multiple runs in multiprocessing
    else:
        print("Invalid output length")

    # Post-process the results
    # Calculate the uncertainty bounds
    glue_experiment.finalize(all_results)
        
    print(f'#### End of the GLUE experiment - Postrun ###')

if __name__ == '__main__':
    main(config_path='3_GLUE_postrun/config.ini')

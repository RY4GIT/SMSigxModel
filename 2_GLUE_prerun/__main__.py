# A main module to run various analysis with CFE model

# Import libraries
import os
import multiprocessing as mp
import configparser
from agent import Agent_GLUE_CFE

def read_config(config_path): 
    config = configparser.ConfigParser()
    config.read(config_path)
    return config

def main(config_path):

    print(f'#### Start GLUE experiment - Prerun ###')
    
    # Read config
    config = read_config(config_path=config_path)
    
    # Create a GLUE instance
    glue_experiment = Agent_GLUE_CFE(config)

    # Sample parameters
    sampled_param_sets = glue_experiment.generate_params()
    
    # Just to run one model run for a test, 
    if config['Multiprocessing']['run_multithreadding'] == 'False':
        all_results = glue_experiment.simulation(sampled_param_set=sampled_param_sets[0])
    
    # Start multiple runs in multiprocessing
    elif config['Multiprocessing']['run_multithreadding'] == 'True':
        # Set MP settings
        print("Number of cpu : ", mp.cpu_count())
        pool_nprocess = int(config['Multiprocessing']['pool_nprocess'])
        os.environ['NUMEXPR_MAX_THREADS'] = config['Multiprocessing']['NUMEXPR_MAX_THREADS']
        
        # Run MP 
        pool = mp.Pool(processes=pool_nprocess)
        all_results = pool.map(glue_experiment.simulation, sampled_param_sets)
        pool.close()
        pool.join()
        print(f"Finished GLUE runs (n = {config['GLUE']['nrun']})")

    # Post-process the results
    glue_experiment.finalize(all_results)

if __name__ == '__main__':
    main(config_path='2_GLUE_prerun/config.ini')
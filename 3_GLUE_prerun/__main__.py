# A main module to run various analysis with CFE model

# Import libraries
import os
import multiprocessing as mp
import configparser
from agent import Agent_GLUE_CFE
import time


def read_config(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def main(config_path):
    start_time = time.time()
    print(f"#### Start GLUE experiment - Prerun ###")

    # Read config
    config = read_config(config_path=config_path)

    # Create a GLUE-CFE agent
    glue_experiment = Agent_GLUE_CFE(config)

    # Sample random parameters
    sampled_param_sets = glue_experiment.generate_params()

    # Just to excute a single model run for a test
    if config["Multiprocessing"]["run_multithreadding"] == "False":
        print("--- Single run mode ---")
        all_results = [
            glue_experiment.simulation(sampled_param_set=sampled_param_sets[0])
        ]

    # To execute multiple runs in multiprocessing
    elif config["Multiprocessing"]["run_multithreadding"] == "True":
        print("--- Multithreadding mode ---")

        # MP settings
        print("Number of cpu : ", mp.cpu_count())
        pool_nprocess = int(config["Multiprocessing"]["pool_nprocess"])
        os.environ["NUMEXPR_MAX_THREADS"] = config["Multiprocessing"][
            "NUMEXPR_MAX_THREADS"
        ]

        # Run MP
        print("--- Simulation started ---")
        pool = mp.Pool(processes=pool_nprocess)
        all_results = pool.map(glue_experiment.simulation, sampled_param_sets)
        pool.close()
        pool.join()
        print(f"--- Finished GLUE runs (n = {config['GLUE']['nrun']}) ---")

    # Post-process the results
    glue_experiment.finalize(all_results=all_results)

    print(f"#### End of the GLUE experiment - Prerun ###")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")


if __name__ == "__main__":
    main(config_path="3_GLUE_prerun/config.ini")

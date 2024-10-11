# A main module to run various analysis with CFE model
# To implement sensitivity analysis with SALib. Currently this module supports Morris and Sobol analysis

# Import libraries
import multiprocessing as mp
import os
from agent import Agent_SALib_CFE
import configparser
import time
from tqdm import tqdm


def map_with_progress(func, iterable, num_processes):
    with mp.Pool(num_processes) as pool:
        results = list(tqdm(pool.imap(func, iterable), total=len(iterable)))
    pool.close()
    pool.join()
    return results


def main():
    start = time.perf_counter()

    config = configparser.ConfigParser()
    config.read("2_sensitivity_analysis/config.ini")

    print(f"### Start {config['SALib']['method']} sensitivity analysis ###")

    # Prepare temporary file location
    source_path = config["PATHS"]["cfe_config"]
    directory = os.path.join(
        os.path.dirname(source_path),
        "temporary_parameter_files_for_sensitivity_analysis",
    )
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Implementation
    salib_experiment = Agent_SALib_CFE(config=config)
    sampled_param_sets = salib_experiment.sample()

    # Start of the multiprocessing
    # MP settings
    print("Number of cpu : ", mp.cpu_count())
    pool_nprocess = int(config["Multiprocessing"]["pool_nprocess"])
    os.environ["NUMEXPR_MAX_THREADS"] = config["Multiprocessing"]["NUMEXPR_MAX_THREADS"]

    # Run MP
    print("--- Evaluation started ---")
    # pool = mp.Pool(processes=pool_nprocess)
    results = map_with_progress(
        salib_experiment.simulation, sampled_param_sets, pool_nprocess
    )
    # results = pool.map(salib_experiment.simulation, sampled_param_sets)
    print(f"--- Finished evaluation runs ---")

    salib_experiment.finalize(results)

    end = time.perf_counter()
    print(f"Run took : {(end - start):.6f} seconds")


if __name__ == "__main__":
    main()

# A main module to run various analysis with CFE model
# To implement sensitivity analysis with SALib. Currently this module supports Morris and Sobol analysis

# Import libraries
import sys
# if not sys.warnoptions:
#     import warnings
#     warnings.simplefilter("ignore")
from agent import Agent_SALib_CFE
import configparser
import time

def main():
    
    start = time.perf_counter()
    
    config = configparser.ConfigParser()
    config.read('1_sensitivity_analysis/config.ini')
    
    print(f"### Start {config['SALib']['method']} sensitivity analysis ###")

    # Implementation
    salib_experiment = Agent_SALib_CFE(config=config)
    salib_experiment.run()
    salib_experiment.finalize()
    
    end = time.perf_counter()
    print(f"Run took : {(end - start):.6f} seconds")
        
if __name__ == '__main__':
    main()

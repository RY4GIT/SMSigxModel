import os
import sys
sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs")
import cfe
import spotpy_cfe
import spotpy
import numpy as np

# specify current directory create output directory if it does not exist
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model")
os.getcwd()
out_file_path = '../4_out/test/'
if not os.path.exists(out_file_path):
    os.mkdir(out_file_path)
data_file_path = '../2_data_input/test'

def main():

    # ensambles ==================================================
    cfe1 = cfe.CFE(os.path.join(data_file_path, 'cat_58_config_cfe.json'))
    # cfe1.run_unit_test()

    # sensitivity analysis ==================================================
    # Initialize
    spot_setup = spotpy_cfe.spot_setup(cfe=cfe1)
    out_fn_sa = out_file_path + 'results'

    # Select number of maximum repetitions
    # CHeck out https://spotpy.readthedocs.io/en/latest/Sensitivity_analysis_with_FAST/
    # How to determine an appropriate number of repetitions
    rep = 5

    # Start a sensitivity analysis
    sampler = spotpy.algorithms.fast(spot_setup, dbname=out_fn_sa, dbformat='csv', db_precision=np.float32, save_sim = False)
    sampler.sample(rep)

    # Load the results gained with the fast sampler
    results = spotpy.analyser.load_csv_results(out_fn_sa)

    # Example plot to show the sensitivity index of each parameter
    spotpy.analyser.plot_fast_sensitivity(results, number_of_sensitiv_pars=3)

    # Example to get the sensitivity index of each parameter
    SI = spotpy.analyser.get_sensitivity_of_fast(results)

    """
    cfe1.initialize()
    cfe1.update()
    cfe1.update_until(4)
    cfe1.finalize()
    """

if __name__ == '__main__':
    main()

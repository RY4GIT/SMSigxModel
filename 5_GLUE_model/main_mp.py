# A main module to run various analysis with CFE model

# Import libraries
import os
import sys
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/5_GLUE_model")

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
from glue_cfe_mp import MyGLUE

def main(out_path='', config_path_CFE='', config_path_GLUE='', nrun=1, eval_criteria=dict()):
    # To implement sensitivity analysis with GLUE.
    print('Start GLUE analysis')

    # Create output directry if not exist
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # Create a GLUE instance
    glue_instance = MyGLUE(
        out_path=out_path,
        config_path_CFE=config_path_CFE,
        nrun=nrun
    )

    # Read parameter bounds defined from an excel file
    df = pd.read_excel(config_path_GLUE)
    df_param_to_calibrate = df[df['calibrate'] == 1]
    params = len(df_param_to_calibrate) * [None]
    for i in range(len(params)):
        params[i] = spotpy.parameter.Uniform(
            df_param_to_calibrate['name'][i],
            low=df_param_to_calibrate['lower_bound'][i],
            high=df_param_to_calibrate['upper_bound'][i]
        )

    # Generate random parameters
    sampled_params = nrun * [None]
    for i in range(len(sampled_params)):
        sampled_params[i] = [i, spotpy.parameter.generate(params)]

    # glue_instance.simulation(sampled_params=sampled_params[0])
    # Start multiple runs in multiprocessing
    pool = mp.Pool(processes=3)
    all_results = pool.map(glue_instance.simulation, sampled_params)
    pool.close()
    pool.join()

    # Post-process the results
    glue_instance.save_results_to_df(all_results)

    # Output the results
    glue_instance.to_csv(df_param_to_calibrate=df_param_to_calibrate)

    print(f'Finished GLUE run ({nrun} runs)')
    print(f'Saved results to: {out_path}')

if __name__ == '__main__':

    # measure the running time
    measuretime = False
    if measuretime:
        pr = cProfile.Profile()
        pr.enable()

    # ===============================================
    # =========== GLUE ANALYSIS ==============
    # ===============================================


    # The main run
    main(
        out_path='../6_out/Mahurangi/ex111',
        config_path_CFE='../2_data_input/Mahurangi/parameters/config_cfe_template.json',
        config_path_GLUE='../2_data_input/Mahurangi/parameters/ex1_GLUE_config.xlsx',
        nrun=3
    )

    # measure the time
    if measuretime:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        # print(s.getvalue())
        ps.dump_stats('runtime.txt')

    # "C:\Users\flipl\anaconda3\envs\CFE\Scripts\snakeviz.exe"
    # snakeviz "G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\3_model\runtime.txt"
    # visualize to type the above in terminal
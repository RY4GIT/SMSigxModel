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
from glue_post_analysis import MyGLUEPost

def main(out_path='', config_path_CFE='', path_GLUE_output=''):
    # To reproduce behavioral runs and evaluate
    print("Post evaluation of GLUE analysis")

    # Create output directry if not exist
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # Parameter bounds defined from an excel file
    df_params = pd.read_csv(os.path.join(path_GLUE_output, 'paramter_priori.xls'), index_col=0)
    df_glue_results = pd.read_csv(os.path.join(path_GLUE_output, 'glue_results.xls'), index_col=0)
    behavioral_idx = df_glue_results.index[df_glue_results['Behavioral'].values]
    df_params_behavioral = df_params.iloc[behavioral_idx].copy()
    print(f"Total {len(behavioral_idx)} runs")

    # Create a GLUE instance
    glue_instance = MyGLUEPost(
        out_path=out_path,
        config_path_CFE=config_path_CFE
    )

    # Generate random parameters
    behavioral_params = len(behavioral_idx) * [None]
    for i, run_id in enumerate(behavioral_idx):
        behavioral_params[i] = [run_id, df_params_behavioral.loc[run_id]]

    # all_results = glue_instance.simulation(behavioral_param=behavioral_params[0])

    # Start multiple runs in multiprocessing
    pool = mp.Pool(processes=4)
    all_results = pool.map(glue_instance.simulation, behavioral_params)
    # all_results = pool.map(glue_instance.simulation, behavioral_params[0:3])
    pool.close()
    pool.join()

    # Post-process the results
    glue_instance.save_results_to_df(all_results)

    # Output the results
    glue_instance.to_csv()

    print(f'Finished the run')
    print(f'Saved results to: {out_path}')

if __name__ == '__main__':

    main(
        out_path= r'../8_out/Mahurangi/ws2_ex2',
        config_path_CFE= r'../2_data_input/Mahurangi/parameters/config_cfe_template.json',
        path_GLUE_output= r'../6_out/Mahurangi/ws2_ex2'
    )

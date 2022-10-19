# A main module to run various analysis with CFE model

# Import libraries
import os
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
import cProfile, pstats, io
from pstats import SortKey
import multiprocessing as mp
import snakeviz


# https://docs.python.org/3/library/profile.html#module-cProfile

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/libs/")
from glue_cfe import MyGLUE

# Specify current directory create output directory if it does not exist
os.chdir("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/5_GLUE_model")
os.getcwd()

def main(out_path='', config_path_CFE='', config_path_GLUE='', nrun=1, eval_criteria=dict()):
    # To implement sensitivity analysis with GLUE.
    print('Start GLUE analysis')

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # Start GLUE
    glue_instance = MyGLUE(
        out_path=out_path,
        config_path_CFE=config_path_CFE,
        config_path=config_path_GLUE,
        nrun=nrun,
        eval_criteria=eval_criteria
    )

    glue_instance.simulation()
    glue_instance.post_process()

    # Output the results
    glue_instance.to_csv()
    glue_instance.plot(plot_type="dotty")
    glue_instance.plot(plot_type="dotty_interaction")
    glue_instance.plot(plot_type="param_hist")
    glue_instance.plot(plot_type="timeseries")

if __name__ == '__main__':

    # measure the running time
    measuretime = True
    if measuretime:
        pr = cProfile.Profile()
        pr.enable()


    # ===============================================
    # =========== GLUE ANALYSIS ==============
    # ===============================================
    eval_criteria_NSE = {
        0: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'NSE', 'threshold': 0.5}
    }

    eval_criteria_KGE = {
        0: {'variable_to_analyze': 'Flow', 'metric': 'KGE', 'threshold': 0.5}
    }

    eval_criteria_season = {
        0: {'variable_to_analyze': 'Soil Moisture Content', 'metric': 'season_transition', 'threshold': 30}
    }


    # variable_to_analyze: ["Flow", "Soil Moisture Content"]
    # metric = ["NSE", "KGE", "season_transition"]

    main(
        out_path='../6_out/Mahurangi/ex3',
        config_path_CFE='../2_data_input/Mahurangi/parameters/ex1_config_cfe.json',
        config_path_GLUE='../2_data_input/Mahurangi/parameters/ex1_GLUE_config.xlsx',
        nrun=200,
        eval_criteria=eval_criteria_season
    )

    # main(
    #     out_path='../6_out/Mahurangi/ex3',
    #     config_path_CFE='../2_data_input/Mahurangi/parameters/ex1_config_cfe.json',
    #     config_path_GLUE='../2_data_input/Mahurangi/parameters/ex1_GLUE_config.xlsx',
    #     nrun=200,
    #     eval_criteria=eval_criteria_seasonsig
    # )


    # process1 = mp.Process(target=main, kwargs={'out_path': '../6_out/Mahurangi/ex1',
    #                                             'config_path_CFE': '../2_data_input/Mahurangi/parameters/ex1_config_cfe.json',
    #                                              'config_path_GLUE': '../2_data_input/Mahurangi/parameters/ex1_GLUE_config.xlsx',
    #                                              'nrun':3,
    #                                               'eval_criteria':eval_criteria})
    #
    # process2 = mp.Process(target=main, kwargs={'out_path': '../6_out/Mahurangi/ex1',
    #                                            'config_path_CFE': '../2_data_input/Mahurangi/parameters/ex1_config_cfe.json',
    #                                            'config_path_GLUE': '../2_data_input/Mahurangi/parameters/ex1_GLUE_config.xlsx',
    #                                            'nrun': 3,
    #                                            'eval_criteria': eval_criteria})
    #
    # process3 = mp.Process(target=main, kwargs={'out_path': '../6_out/Mahurangi/ex1',
    #                                             'config_path_CFE': '../2_data_input/Mahurangi/parameters/ex1_config_cfe.json',
    #                                              'config_path_GLUE': '../2_data_input/Mahurangi/parameters/ex1_GLUE_config.xlsx',
    #                                              'nrun':3,
    #                                               'eval_criteria':eval_criteria})
    #
    # process4 = mp.Process(target=main, kwargs={'out_path': '../6_out/Mahurangi/ex1',
    #                                             'config_path_CFE': '../2_data_input/Mahurangi/parameters/ex1_config_cfe.json',
    #                                              'config_path_GLUE': '../2_data_input/Mahurangi/parameters/ex1_GLUE_config.xlsx',
    #                                              'nrun':3,
    #                                               'eval_criteria':eval_criteria})
    # process2 = mp.Process(target=main)

    # process1.start()
    # process2.start()
    # process3.start()
    # process4.start()


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
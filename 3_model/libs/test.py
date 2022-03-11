import pandas as pd
compare_results_file = "G:\\Shared drives\\Ryoko and Hilary\\SMSigxModel\\analysis\\2_data_input\\Mahurangi\\short\\test_sm_basinavg.csv"
unit_test_data = pd.read_csv(cfe1.compare_results_file)
data = pd.DataFrame().reindex_like(unit_test_data)
import time
import numpy as np
import pandas as pd
import json
import sys

sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/libs/cfe/py_cfe")
import bmi_cfe

cfe = bmi_cfe.BMI_CFE('G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/2_data_input/unit_test/config_cfe.json')
cfe.initialize()
cfe.run_unit_test(print_fluxes=True, warm_up=True)
print(cfe.cfe_output_data)
cfe.finalize(print_mass_balance=True)

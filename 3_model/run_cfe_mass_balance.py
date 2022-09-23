import time
import numpy as np
import pandas as pd
import json
import sys
sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs/cfe/py_cfe")
import cfe
import bmi_cfe

cfe1 = bmi_cfe.BMI_CFE('G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/2_data_input/debug/config_cfe.json')
cfe1.initialize()
cfe1.run_unit_test(print_fluxes=True)
print(cfe1.cfe_output_data)
cfe1.finalize(print_mass_balance=True)

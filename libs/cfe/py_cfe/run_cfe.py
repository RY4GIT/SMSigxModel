import time
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import bmi_cfe

""" 
    initialize: Perform startup tasks for the model.
    update: Advance model state by one time step. Calls the function run_cfe from cfe.py
    update_until: Advance model state until the given time.
    finalize: Perform tear-down tasks for the model.
    get_value: Get a copy of values of a given variable.
    set_value: Set the values of a given variable.
    etc.
    These functions need to be called by a framework or driving code, an example of which is below.
"""
"""
# Create an instance of the model with a specific configuration that corresponds to a particular catchmenmt.
cfe_instance = bmi_cfe.BMI_CFE('./cat_58_config_cfe.json')

# This initialization function should perform all tasks that are to take place before entering the model’s time loop. Models should be refactored, if necessary, to read their inputs (which could include filenames for other input files) from a configuration file. BMI does not impose any constraint on how configuration files are formatted.
cfe_instance.initialize()

# Open the forcing file contained within the configuration file. We can run the model with any forcing. This is only an example. The path to the forcing file is contained within the configuration file, but it doesn't really need to be. This is just for organization.
with open(cfe_instance.forcing_file, 'r') as f:
    df_forcing = pd.read_csv(f)

# We will want to visualize the model output
outputs=cfe_instance.get_output_var_names()
output_lists = {output:[] for output in outputs}

# Now we loop through the forcing data and use it to run the model at each time step
for precip in df_forcing['APCP_surface']:

    cfe_instance.set_value('atmosphere_water__time_integral_of_precipitation_mass_flux', precip)

    cfe_instance.update()

    for output in outputs:
        output_lists[output].append(cfe_instance.get_value(output))

# The finalize function should perform all tasks that take place after exiting the model’s time loop. This typically includes deallocating memory, closing files and printing reports.
cfe_instance.finalize(print_mass_balance=True)

istart_plot=490
iend_plot=550
x = list(range(istart_plot, iend_plot))
for output in outputs:
    plt.plot(x, output_lists[output][istart_plot:iend_plot], label=output)
    plt.title(output)
    plt.legend()
    plt.show()
    plt.close()
"""
# Here we are just going to run the unit test that compares with the origional author code. Kind of the same thing done above, but there is a function in the BMI code that does it all at once.
cfe_instance2 = bmi_cfe.BMI_CFE('./cat_58_config_cfe.json')
cfe_instance2.initialize()
cfe_instance2.run_unit_test(plot=True, print_fluxes=False)
cfe_instance2.finalize()
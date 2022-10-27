import shutil
import os
n=1

original = r"..\2_data_input\Mahurangi\parameters\config_cfe_template.json"
target = os.path.join(r"..\2_data_input\Mahurangi\parameters", f"config_cfe_template_{n}.json")

shutil.copyfile(original, target)
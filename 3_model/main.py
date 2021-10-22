import os
import sys
sys.path.append("G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model/libs")
import cfe

# specify current directory create output directory if it does not exist
os.chdir('G://Shared drives/Ryoko and Hilary/SMSigxModel/analysis/3_model')
os.getcwd()
out_file_path = '../4_out/test'
if not os.path.exists(out_file_path):
    os.mkdir(out_file_path)

data_file_path = '../2_data_input/test'

def main():

    # ensambles
    cfe1 = cfe.CFE(os.path.join(data_file_path, 'cat_58_config_cfe.json'))

    cfe1.run_unit_test()
    """
    cfe1.initialize()
    cfe1.update()
    cfe1.update_until(4)
    cfe1.finalize()
    """

    """
    original ... 
    cfe1 = cfe.CFE(os.path.join(data_file_path, 'cat_58_config_cfe.json'))
    cfe1.run_unit_test()
    cfe1.initialize()
    cfe1.update()
    cfe1.update_until(4)
    cfe1.finalize()
    """
if __name__ == '__main__':
    main()

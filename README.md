# SMSigxModel
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![versions](https://img.shields.io/pypi/pyversions/hydra-core.svg) [![CodeStyle](https://img.shields.io/badge/code%20style-Black-black)]()  
This project explores the application of soil moisture signature ([Branger et al., 2019](https://doi.org/10.1002/hyp.13645); [Araki et al., 2020](https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.14553)) to enhance streamflow and soil moisture prediction within a rainfall-runoff model.

This code is currently in development and does not guarantee to run yet. 

## Installation
Use conda to create your own env based on our ```environment.yml``` file

```bash
conda env create -f environment.yml
conda activate CFE
```

## Running this code
Each directory houses a package designed to run each experiment step. These packages utilize the `__main__.py` functions, and they come with configuration files (typically named `config.ini`) along with some additional configuration files.

The primary branch is set up to run the Mahurangi test case. If you wish to use a different case, adjust the configuration files accordingly.

### 0_data_preprocessing
Formats Mahurangi & Little Washita data. Formatted datasets are saved in the `data` folder. If you only plan to run experiments #1 through #4, this step isn't required.

### 1_sensitivity_analysis
Conducts a Morris sensitivity analysis.

To run:
```bash
cd 1_sensitivity_analysis
python __main__.py
```

- ```config.ini```
    - The main config file. Copy `example_config.ini` and change it to your desired setting
- ```config_SALib.csv```
    - The configuration file that defines the calibrated parameter and its bounds

### 2_GLUE_prerun
Runs the rainfall-runoff model with various randomly-generated parameters in preparation for the next step. This step is not dependent on the `1_sensitivity_analysis`; If you are only interested in the GLUE experiments, start from this step.  

To run:
```bash
cd 2_GLUE_prerun
python __main__.py
```

- ```config.ini```
    - The main config file. Copy `example_config_xxx.ini` and change it to your desired setting
- ```config_GLUE.csv```
    - The configuration file that defines the calibrated parameter and its bounds

### 3_GLUE_postrun 
Analyzes the output from `2_GLUE_prerun` using GLUE.

To run:
```bash
cd 3_GLUE_postrun
python __main__.py
```

- ```config.ini```
    - The main config file. Copy `example_config_xxx.ini` and change it to your desired setting
- ```.3_GLUE_postrun\config_criteria\GLUE_criteria_{criteria_id}.json```
    - The configuration file that defines the criteria of GLUE behavioral threshold. The criteria_id in the filename corresponding to config["GLUE"]["criteria_id"]

### 4_post_analysis
Executes post-analysis and visualizes the results from 3_GLUE_postrun. Follow the Jupyter Notebooks in numerical order for this step.

## Data
- Input data requirements
  - `time`: hourly time interval
  - `precip_rate`: preipictation in [m/hr]
  - `PET`: potential evapotranspiration in [m/hr]
- Observation data requirements
  - `Time`: hourly time interval
  - `Rainfall`: preipictation in [m/hr]
  - `Flow`: streamflow discharge in [m/hr]
  - `Soil Moisture Content`: volumetrics soil water content in [m3/m3]
  - `Direct runoff` `GIUH Runoff` `Lateral Flow` `Base Flow`: not used in analysis but BMI_CFE requires it. All values set `0` as default. 
  - `Total Discharge`: streamflow discharge in [m3/s] (`Flow` multiplied by area and time conversion)

## Resources

#### Sensitivity analysis 
SALib documentation: https://salib.readthedocs.io/en/latest/

    Iwanaga, T., Usher, W., & Herman, J. (2022). Toward SALib 2.0: Advancing the accessibility and interpretability of global sensitivity analyses. Socio-Environmental Systems Modelling, 4, 18155. doi:10.18174/sesmo.18155

    Herman, J. and Usher, W. (2017) SALib: An open-source Python library for sensitivity analysis. Journal of Open Source Software, 2(9). doi:10.21105/joss.00097

#### GLUE analysis
    Beven KJ, Binley AM. 1992. The future of distributed models: model calibration and uncertainty prediction. Hydrological Processes 6: 279–298. 

#### A rainfall-runoff model
The CFE model is designed to be a simplified and functionaly equivalent model of the National Water Model. The model code was originally written by Dr. Fred Ogden and converted to BMI-compliant format in the Next-Gen framework by NOAA-OWP. The official CFE code by Dr. Fred Oden and NOAA-OWP lives [here](https://github.com/NOAA-OWP/cfe/).  [The Python version of the code](https://github.com/NWC-CUAHSI-Summer-Institute/cfe_py) is developed for the prototyping, research, and development. This code is developed upon the Python version and for research purpose. 

### License
[MIT](https://choosealicense.com/licenses/mit/)
Project Template provided by: https://github.com/moemen95/Pytorch-Project-Template

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Citation
To cite this code, email the author raraki8159 at sdsu dot edu
# SMSigxModel
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![versions](https://img.shields.io/pypi/pyversions/hydra-core.svg) [![CodeStyle](https://img.shields.io/badge/code%20style-Black-black)]()
To be edited.

#### Conceptual Functional Equivalent (CFE) Model
The CFE model is designed to be a simplified and functionaly equivalent model of the National Water Model. The model code was originally written by Dr. Fred Ogden and converted to BMI-compliant format in the Next-Gen framework by NOAA-OWP. The official CFE code by Dr. Fred Oden and NOAA-OWP lives [here](https://github.com/NOAA-OWP/cfe/).  [The Python version of the code](https://github.com/NWC-CUAHSI-Summer-Institute/cfe_py) is developed for the prototyping, research, and development. This code is developed upon the Python version and for research purpose. 

## Installation
Use conda to create your own env based on our ```environment.yml``` file

```bash
conda env create -f environment.yml
conda activate CFE
```

## Running this code
To be edited.

The main branch code is currently configured to run the xxxx test case. If you want to use your own case, you will need to manage three config files located here:

- ```xxx```
    - The main config file. Choose appropriate xxxx

To run the code, just run the following command inside the xxxx folder:

```python  .```

## Reference
SALib documentation: https://salib.readthedocs.io/en/latest/
    Iwanaga, T., Usher, W., & Herman, J. (2022). Toward SALib 2.0: Advancing the accessibility and interpretability of global sensitivity analyses. Socio-Environmental Systems Modelling, 4, 18155. doi:10.18174/sesmo.18155

    Herman, J. and Usher, W. (2017) SALib: An open-source Python library for sensitivity analysis. Journal of Open Source Software, 2(9). doi:10.21105/joss.00097


## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)

Project Template provided by: https://github.com/moemen95/Pytorch-Project-Template
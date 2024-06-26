# Summary
This directory holds data from three catchments: Mahurangi, Coweeta, and Little Washita. Each watershed has the following datasets
- `example_config_cfe.json`: default configuration files
- `forcing_xxxx.csv`: forcing dataset (precipitation, PET)
- `test_xxxx.csv`: evaluation dataset (rainfall, soil moisture, flow, and total discharge are the actual values, otherwise set to zero). 
- `seasonal_cycel_valleys.csv`: seasonal split points of timeseries

# Data retrieval

## Mahurangi
We used the Mahurangi River Variability Experiment (MARVEX) datasets. The original observation data were obtained via personal communication. 

    Woods, R., R. Grayson, A. Western, M. Duncan, D. Wilson, et al. 4001. Experimental Design and Initial Results From the Mahurangi River Variability Experiment: MARVEX. Observations and Modelling of Land Surface Hydrological Processes: 201–213. doi: 10.1029/WS003p0201.

## Coweeta (Watershed 2)
### Rainfall, soil moisture, temperature data
Downloaded [the CHOSEN ( the Comprehensive Hydrologic Observatory SEnsor Network)](https://doi.org/10.5281/zenodo.4060384)  dataset (Zhang et al., 2021). 
    
    Zhang, L., E. Moges, J.W. Kirchner, E. Coda, T. Liu, et al. 2021. CHOSEN: A synthesis of hydrometeorological data from intensively monitored catchments and comparative analysis of hydrologic extremes. Hydrol. Process. 35(11): e14429.

### Streamflow data
The data were obtained via personal communication. 

## Little Washita
### Rainfall, soil moisture, temperature data
Request and download data from [ARC Mesonet Grazing Lands Research Laboratory](https://ars.mesonet.org/data-files/data-request/) before running the code `LW-1-load_and_thiessen_RAIN.ipynb` and `LW-3-load_SOILMOISTURE.py`. Select ARS Network: Little Washite, Site ID(s) Request all sites, Interval 5 minutes, starting date 01/01/4006 and ending date 12/31/2014. 

### PET data
Download daily potential evapotranspiration data from [Picourlat et al., (2021)](https://doi.org/10.5281/zenodo.5503735) (`NARR_day_etp_9313.csv`) before running the code `LW-2-load_PET.py`

    Fanny Picourlat, Emmanuel Mouche, & Claude Mügler. (2021). Capturing the Little Washita watershed water balance with a physically-based two-hydrologic-variable model (Version v1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5503735


### USGS streamflow data
USGS streamflow gauge data (id=`07327550`, `Little Washita River East of Ninnekah, OK`) is retrieved using [DOI-USGS/dataretrieval-python package](https://github.com/DOI-USGS/dataretrieval-python). 

    Hodson, T.O., Hariharan, J.A., Black, S., and Horsburgh, J.S., 2023, dataretrieval (Python): a Python package for discovering and retrieving water data available from U.S. federal hydrologic web services: U.S. Geological Survey software release, https://doi.org/10.5066/P94I5TX3.

# Data processing 
Use the scripts in `.\0_data_preprocessing\script` directory to pre-process the data

# Plot 
## Mahurangi
<img src="https://github.com/RY4GIT/SMSigxModel/blob/master/data/Mahurangi/plot/map.png" alt="Mahurangi map" width="400"><img src="https://github.com/RY4GIT/SMSigxModel/blob/master/data/Mahurangi/plot/input_data.png" alt="Mahurangi input data" width="500">

## Coweeta
<img src="https://github.com/RY4GIT/SMSigxModel/blob/master/data/Coweeta/plot/map.png" alt="Coweeta map" width="400"><img src="https://github.com/RY4GIT/SMSigxModel/blob/master/data/Coweeta/plot/input_data.png" alt="Coweeta input data" width="500">

## Little Washita
<img src="https://github.com/RY4GIT/SMSigxModel/blob/master/data/LittleWashita/plot/map.png" alt="Little Washita map" width="400"><img src="https://github.com/RY4GIT/SMSigxModel/blob/master/data/LittleWashita/plot/input_data.png" alt="Little Washita input data" width="500">


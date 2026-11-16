# Quantile-Based Diagnostics for Conditioned Citizen Science Phenological Distributions

## Overview

This repository contains the datasets and R scripts supporting the manuscript *“Quantile-Based Diagnostics for Conditioned Citizen Science Phenological Distributions.”*

The example code can be executed in the R environment. Due to the data access policies, the climate data used in this research must be obtained through an application to the Taiwan Climate Change Projection Information and Adaptation Knowledge Platform (TCCIP) (https://tccip.ncdr.nat.gov.tw/index_eng.aspx).

To derive CSCI, please refer to https://github.com/brainbrian2000/Temperature_based_phenological_indicator_under_data_limited_condiction_ERA5_and_citizen_science, which hosts the code and data for the paper titled “Deriving exploratory temperature-based phenological indicators under data-limited conditions: integrating ERA5 and citizen science”.

The repository includes simulated datasets, species-specific climate and phenological observations, and the R code used for QRF modeling and validation.

## Repository Structure

```text
.
├── QRF_phenology_analysis.R
│
├── Simulation/
│   ├── grid_for_QRF_analysis_overlapping_regions.csv
│   └── QRF_simulation_analysis.R
│
└── Data/
    ├── Hyla_chinensis/
    │   ├── VIP1/
    │   │   └── CSCI_Features_for_R_QRF.csv
    │   └── VIP08/
    │       └── CSCI_Features_for_R_QRF.csv
    │
    ├── Zhangixalus_moltrechti/
    │   ├── VIP1/
    │   └── VIP08/
    │
    ├── Fejervarya_limnocharis/
    │   ├── VIP1/
    │   └── VIP08/
    │
    └── Microhyla_fissipes/
        ├── VIP1/
        └── VIP08/
```


## 1. Simulation

The `Simulation/` directory contains a synthetic gridded dataset generated according to the simulation method described in the manuscript.

The simulated grid represents spatial variation in climatic conditions and associated phenological responses. The directory also includes an R script for cross-region QRF training and validation, allowing model performance to be evaluated when training and testing observations originate from different climatic regions.


## 2. Empirical Data

The `Data/` directory contains four species-specific subdirectories:

- *Hyla chinensis*
- *Zhangixalus moltrechti*
- *Fejervarya limnocharis*
- *Microhyla fissipes*

Each species directory contains two subdirectories, `VIP1/` and `VIP08/`, representing different variable selection thresholds based on Variable Importance in Projection (VIP).

Each corresponding `CSCI_Features_for_R_QRF.csv` file contains:

- Phenological observation records used in the analysis.
- Geographic coordinates (latitude and longitude) of the original observation locations.
- Observation dates and day-of-year (DOY) values.
- Calculated CSCI variables used as climatic predictors in QRF modeling.


## 3. Main Analysis Script

The main analysis script, `QRF_phenology_analysis.R`, is located in the repository's root directory.

It contains the analytical workflow for:

- Fitting Quantile Regression Forest models.
- Estimating the 25th, 50th, and 75th conditional quantiles.
- Evaluating in-sample model performance.
- Performing random cross-validation (RCV).
- Performing spatial block cross-validation (SCV).
- Calculating ICP, MPIW, MAE, and pinball loss.


## 4. Reproducibility

The datasets and scripts are organized to facilitate inspection and reproduction of the simulation and empirical QRF analyses.

To reproduce an empirical analysis:

1. Select one of the four species.
2. Select the corresponding `VIP1/` or `VIP08/` dataset.
3. Set the input CSV path in `QRF_phenology_analysis.R`.
4. Install the required R packages, including `quantregForest`, `dplyr`, `tidyr`, `ggplot2`, `sf`, and `blockCV`.
5. Run the script to obtain model predictions and evaluation metrics.

Random seeds are specified in the analysis scripts where applicable to support reproducibility.


## 5. Data Sources and Citation

The empirical analyses use citizen science phenological observations and derived climatic indicators. Details of the data sources, processing procedures, CSCI calculations, and model configurations are provided in the accompanying manuscript.

When using the datasets or code, please cite the associated publication.

## Contact

For questions regarding the datasets or analysis procedures, please contact the corresponding author of the associated manuscript.

# CAFA 5 Kaggle Competition
---

## Overview
This repo is associated with a Kaggle competition. The goal of the competition is to predict the gene ontologies for a set of proteins. 
It contains data ETL, data analysis, feature engineering, and modeling.


## Setup/Installation
1. Presupposes an existing miniconda/conda setup and an environment with `python>3.10` (via `pyproject.toml`)
2. Clone repo and navigate to the root level of the directory housing the repository
3. Within the appropriate conda/python environment, `pip install` or `pip install -e .` if you're actively developing


## Usage

### Data Preparation

1. Get data from [kaggle](https://www.kaggle.com/competitions/cafa-5-protein-function-prediction/data) and house in `./data`
2. Open `notebooks/execute_data_flow.ipynb` and run it in its entirety
3. This will produce two serialized datafiles that have been transformed, featured engineered & merged.
   * Datafiles are outputed to the `data/flow_output/` directory

### Data Analysis
**TODO** See `project_overview.md` for now

### Modeling
**TODO** See `project_overview.md` for now

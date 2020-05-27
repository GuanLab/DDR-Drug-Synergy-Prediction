# Merck DDR Drug Synergy Prediction

This repository contains the code of the final model for DDR drug Synergy Prediction Project Wave 1. 

## Dependency

[Python 3.0](https://www.python.org/download/releases/3.0/)

[LightGBM 2.3.2](https://lightgbm.readthedocs.io/en/latest/index.html)

## Data Resource

To run this repository smoothly, three types of feature data (response, feature and geneset) should be placed under the `/data` directory.  

We also included two extra types features, chemical structure and gene network features. The raw data used to construct network features is avaliable from [MouseNet](http://fntm.princeton.edu). 

## Data Preprocess

QC includes the general analysis for the quality of the above three types of data.

## Split Training and Testing Data

```
python test_split.py
```
this will create test and training set splits from the Drug Synergy dataset.

## Train Models and Create SHAP Analysis Plots

Two models were constructed parallelly to evaluate the drug synergy prediction under different circunstances

* Merck_test_by_cell_line : cross indication models
* Merck_test_by_batch: cross batch models

To run both model, first we run

```
python data_process.py
```

This will construct feature matrix for all the monotherapy features.

Then we run

```
python model_training.py
```

This program will carry out four tasks: 
* Train and validate the the final version of the Drug Synergy prediction model.
* Carry our SHAP analysis on the whole test set
* Carry out SHAP analysis on each of the mode-of-action syngergy pairs
* Print dependency plot of different types of molecular markers of the gene of interest.








# Merck DDR Drug Synergy Prediction
##dependency
```
Python 3
LightGBM
```

This repository contains the code of the final model for DDR drug Synergy Prediction Project Wave 1. 
## Data Resource

To run this mondel smoothly, three types of feature data (response, feature and geneset) should be placed under the /data directory.  

we also included two extra types features, chemical structure and gene network features. the raw data to be used to construct network features is avaliable from MouseNet: http://fntm.princeton.edu 

## data preprocess

QC includes the general analysis for the quality of the above three types of data.

### split training and testing data
```
python test_split.py
```
this will create test and training set splits from the Drug Synergy dataset.

## Model training and create SHAP analysis Plots

Merck_test_by_cell_line : cross indication models
Merck_test_by_batch: cross batch models

both to run both model, first we run
```
python data_process.py
```
This will construct feature matrix for all the monotherapy features

then we run:
```
python model_training.py
```
This program will carry out four tasks: 
1) Train and validate the the final version of the Drug Synergy prediction model.
2) Carry our SHAP analysis on the whole test set
3) Carry out SHAP analysis on each of the mode-of-action syngergy pairs
4) Print dependency plot of different types of molecular markders of the gene of interest.








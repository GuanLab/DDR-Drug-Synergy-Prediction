# Merck DDR Drug Synergy Prediction

This repository contains the code of the final model for DDR drug Synergy Prediction Project Wave 1. 

## Dependency

[Python 3.0](https://www.python.org/download/releases/3.0/)

[LightGBM 2.3.2](https://lightgbm.readthedocs.io/en/latest/index.html)

## Data Resource

To run this repository smoothly, three types of feature data (response, feature and geneset) should be placed under the `/data` directory.  

We also included two extra types features, chemical structure and gene network features. The raw data used to construct network features is avaliable from [MouseNet](http://fntm.princeton.edu). 

## Data Preprocess

* 1. `./QC` includes the general analysis for the quality of the above three types of data.
* 2. Quantile normalization on gene expression levels. Run 
```
cd ./QC
python quantile_expression_level.py
```
This will generate quantile normalized feature sets.

## Split Training and Testing Data

```
python test_split.py
```
this will create test and training set splits from the Drug Synergy dataset.

## Train Final Models and Create SHAP Analysis Plots

Since we designed different cross-validation settings, two sets of models were constructed parallelly to evaluate the drug synergy prediction under different circunstances

* Merck_test_by_cell_line : cross indication models
* Merck_test_by_batch: cross batch models

to start training the above two sets of models, first we need to `cd` to the corresponding directoriy:
```
cd Merck_test_by_cell_line
``` 
or
```
cd Merck_test_by_batch
```

### Data preparation

To run both model, first we run

```
python data_process.py
```
This will generate three files or foldera:

* 1. `cleaned_mono.tsv`: this file contains preprocessed monotherapy features, such as `.identifier_batch`, `.metadata_treatment_remarks`, `.response_aoc` and `.metadata_cancer_type`.

* 2. `feature_header.txt`: this file contains the total features would be used in drug synergy prediction model. This file is especially neccessary for SHAP feature importance analysis.

* 3. `./cell_drug_features`: this folder contains the monotherapy features of monotherapy experiments distinguised by specific drug/cell line combinations. This file is for monotherapy feature construction.
This will construct feature matrix for all the monotherapy features.

### Train Synergy Prediction Models

Our LightGBM incorperate the following four types of features:
* Monotherapy experiment efficacy score (AOC)
* Molecular Markers
* Drug Chemical Structure
* Gene Network

To start model training and validation, we run:

```
python model_training.py
```

This program will carry out four tasks: 
* 1. Train Synergy prediction models to predict aoc score (efficacy) and bliss score (synergy), and validate on the validation sets in cross validation. the accuracy will be printed and saved in `cor_aoc.txt` and `cor_bliss.txt`.
* 2. Carry out SHAP analysis of the drug synergy prediction models.
* 3. For each of the mode-of-action combination, such as synergetic treatment of ATMi and ATRi drugs, we carry out model validation and SHAP feature importance analysis as well.
* 4. Print dependency plot of different types of molecular markers of the gene of interest.

### Downstream SHAP feature importance analysis

This part is crucial for analysing potential molecualr biomarkers. 

``` 
python statistical_analysis.py
``` 
This program will carry out analysis in the following aspects:

* 1. Analyse SHAP results by feature category.
* 2. Analyse ddr gene's feature importance by it's molecular features (cnv, snv, lof and exp)
* 3. Analyse tissue-specific ddr gene features.
* 4. Analyse the interaction between priorized biomarkers (ddr genes/gene clusters/ddr signature)

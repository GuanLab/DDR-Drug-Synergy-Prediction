# DDR Drug Synergy Prediction

This repository contains the code for DDR drug Synergy Prediction.

Our LightGBM incorporate the following four types of features:
* Monotherapy experiment efficacy score (AoC)
* Molecular Markers
* Drug Chemical Structure
* Gene Network

## Dependency

1. model implementation:
* [Python 3.0](https://www.python.org/download/releases/3.0/)
* [LightGBM 2.3.2](https://lightgbm.readthedocs.io/en/latest/index.html)

2. to preprocess chemical structure features, the following packages should have been installed:
* [rdkit](https://www.rdkit.org/docs/Install.html)
* [pubchempy 1.0.4](https://pubchempy.readthedocs.io/en/latest/guide/install.html) 
* [openbabel](https://open-babel.readthedocs.io/en/latest/UseTheLibrary/PythonInstall.html)
* [pybel 0.14.10](https://pypi.org/project/pybel/)

3. to retrieve target gene information: 
* [chembl_webresource_client](https://github.com/chembl/chembl_webresource_client)

* Download LINCS dataset for drug-target information in './QC'
  ```
  wget https://lincs.hms.harvard.edu/db/datasets/20000/results?search=&output_type=.csv
  ```
4. downstream feature analysis
* [shap 0.39.0](https://pypi.org/project/shap/)

## Data Structure

This repository contains the following code, where all the rest of data in this project can be generated from:

- `feature`: for data preparation
	- `data`: deposit the drug response data, geneset, chemical structure, molecular readouts etc.
 	- 1. `responses_training.tsv`: drug responses for prediction, including combination therapy and monotherapy. the later will be used as features.
 	- 2. `genesets.tsv`: geneset annotations from literature
 	- 3. `molecular_features.tsv`: molecular features for all cell lines on gene level, including exp, lof, snv, cnv, lof_pat, coh_pat, ddr
	- 4. `chemical_structures.csv`: chemical structure of small molecule drugs in smile format
  
  To run this repository smoothly, the four datasets above should be placed under the `/data` directory.
 	- `QC`: QC and generate features

- `master:`: the code for model training and evaludation
 	- `master_code`: LightGBM model for response prediction
 	- `ensemble`: stacking generation

- `bash.sh`: main command

## Data Storeage, QC and Preprocess
### QC
1. start QC and generate features:
  ```
  cd ./QC
  python QC.py
  ```
  
  `QC.py`: generate following features:
    
    * geneset_qc: 
        -> `../geneset_features/`
    
    * molecular_qc: 
        -> `../molecular_features/`
        -> `../target_gene/ddr_genes.txt`
    
    * response_qc: 
        -> `../monotherapy_features/`
            -> `monotherapy_features.tsv`
        -> `./reproducibility`
        -> `../target_gene/`
            -> `all_drugs_summary.csv`
    
    * encode_categorical_feature:
        -> `../categorical_embedding_features/`
            -> `batch.json`
            -> `cell_line.json`
            -> `moa.json` : mode of action encoding
            -> `treatment.json` : drug name encoding
            -> `remark.json`
            -> `cancer_type.json`
            -> `cancer_subtype.json`

    * cancer_dependency_qc:
        -> `./cancer_dependency_features/`: depMap features; may not use. TODO

2. pull drug targets: 
  ```
  python pull_drug_target.py
  ```
  
  `pull_drug_target.py`: pull drug target information from public databases.
  update `./target_gene/`:
   - 1. the drug target information will be put in `all_drugs_summary.csv`
   - 2. generate a json file `drug2gene.json` for faster target access.

3. split training and testin dataset by cell line/indication/mode-of-action etc.
  ```
  python split_train_test.py 
  ```
  see `split_train_test.py` for details.


## Train Final Models and Create SHAP Analysis Plots

Since we designed different cross-validation settings, two sets of models were constructed parallelly to evaluate the drug synergy prediction under different circunstances.

* Merck_test_by_cell_line : cross indication models
* Merck_test_by_batch: cross batch models

to run the models, run:

```
sh bash.sh
```

for details of training/testing, check the main code by:

```
cd master_code
python main.py -h
```
it would show the usage of main pipeline code:
```
usage: main.py [-h] [-p PATH] [--exclude_synergy_batch] [--exclude_cell_line]
               [--exclude_cancer_type] [-s PREDICTED_SCORE]
               [-f FEATURES [FEATURES ...]]

Build drug synergy prediction machine learning models.

optional arguments:
optional arguments:
  -h, --help            show this help message and exit
  -p PATH, --path PATH  Path to cross-validation split data
  --exclude_synergy_batch
                        If specified, exclude synergy batch as feature. If path is 'test_by_batch', exclude synergy batch automatically.
  --exclude_cell_line   If specified, exclude cell line as feature. If path is 'test_by_cell_line', exclude cell line automatically.
  --exclude_cancer_type
                        If specified, exclude cancer type/subtype as feature. If path is 'test_by_indication', exclude cancer type automatically.
  -s PREDICTED_SCORE, --predicted_score PREDICTED_SCORE
                        Synergy measurement as training target.
                                    aoc: drug efficacy;
                                    bliss: drug synergy;
                                    (default: aoc)
                                    
  -f FEATURES [FEATURES ...], --features FEATURES [FEATURES ...]
                        Feature options in building feature set:
                                    monotherapy: monotherapy responses;
                                    moa: mode-of-action categorical feature;
                                    drug_name: drug name as categorical feature;
                                    molecular: four types of single gene molecular biomarkers,two types of gene cluster biomarkers,78 ddr markers;
                                    target_gene: target gene information (target gene information can only be used when molecular information is included);
                                    network: gene bayesian functional network information.
                                    geneset: 9,450 geneset features;
                                    chemical_structure: drug chemical structure fingerprints;
                                    dependency: DepMap cancer cell line gene dependency features;
                                    synthetic_lethality: synthetic lethality information
                                    (default: [monotherapy, moa, drug_name, molecular, target_gene, geneset, chemical_structure, dependency, synthetic_lethality])
                                    
  --surrogate_gene SURROGATE_GENE
                         number of top genes for surrogate model:
                                    chose one from: 2000, 1000, 500, 250, 200, 125, 100, 75, 50, 25, 20, 10,
                                    if not specified, default as None
                                    
  --surrogate_geneset SURROGATE_GENESET
                         number of top genes for surrogate model:
                                    chose one from: 2000, 1000, 500, 250, 200, 125, 100, 75, 50, 25, 20, 10,
                                    if not specified, default as None
                                    
  --surrogate_chem SURROGATE_CHEM
                         type of chemical stucture feature fingerprint for surroagte model:
                                    chose one from: MACCS, Morgan, FP2, FP3, FP4, RDK
                                    if not specified, default as None
                                    
  --held_out_test       if used, test on the held-out set
  --evaluate_shap EVALUATE_SHAP [EVALUATE_SHAP ...]
                        subset type for subset specific SHAP analysis.
                                    example:
                                    'Train'
                                    'Test'
                                    'moa': Test split by ATM/ATR/DNAPK
                                    'all_moa': Test split by moa pairs
                                    'cell_line'
                                    'tissue'
                                    
  --skip_training       if used, skip training.

```

This program will carry out four tasks: 
* 1. Train Synergy prediction models to predict aoc score (efficacy) and bliss score (synergy), and validate on the validation sets in cross validation. the accuracy will be printed and saved in `cor_aoc.txt` and `cor_bliss.txt`.
* 2. Carry out SHAP analysis of the drug synergy prediction models.
* 3. For each of the mode-of-action combination, such as synergetic treatment of ATMi and ATRi drugs, we carry out model validation and SHAP feature importance analysis as well.
* 4. Print dependency plot of different types of molecular markers of the gene of interest.

### Downstream SHAP feature importance analysis

This part is crucial for analysing potential molecualr biomarkers. 

```
cd master_code
cd downstream_analysis 
python statistical_analysis.py
``` 
This program will carry out analysis in the following aspects:

* 1. Analyse SHAP results by feature category.
* 2. Analyse ddr gene's feature importance by it's molecular features (cnv, snv, lof and exp)
* 3. Analyse tissue-specific ddr gene features.
* 4. Analyse the interaction between priorized biomarkers (ddr genes/gene clusters/ddr signature)

get top genes/genesets for surrgate models:
```
python select_surrogate_features.py
```

calculate improved efficacy against baseline:

```
python calculate_efficacy.py
```

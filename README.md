# Merck DDR Drug Synergy Prediction

This repository contains the code for DDR drug Synergy Prediction.

![Figure1. Overview of the DDR comination drug response model.]('Figure1_model.png')

Our LightGBM incorporate the following four types of features:

* basic drug information: drug name and mode of action
* monotherapy experiment efficacy score (AoC)
* drug chemical structure
* molecular signatures: cnv, snv, exp, lof signatures of ~2,700 DDR genes, coh/lof patterns of gene sets, and 78 ddr signatures
* drug-target interactions
* tissue-specific gene networks
* synthetic lethality

## Dependencies

For building machine learning model:
* [Python 3.0](https://www.python.org/download/releases/3.0/)
* [LightGBM 2.3.2](https://lightgbm.readthedocs.io/en/latest/index.html)

To preprocess chemical structure features, the following packages should have been installed:
* [rdkit](https://www.rdkit.org/docs/Install.html)
* [pubchempy 1.0.4](https://pubchempy.readthedocs.io/en/latest/guide/install.html) 
* [openbabel](https://open-babel.readthedocs.io/en/latest/UseTheLibrary/PythonInstall.html)
* [pybel 0.14.10](https://pypi.org/project/pybel/)

To retrieve target gene information:
* from ChemBL: [chembl_webresource_client](https://github.com/chembl/chembl_webresource_client)
* from LINCS: https://lincs.hms.harvard.edu/db/datasets/20000/results?search=&output_type=.csv
* from DGI: https://dgidb.org/downloads

To retrieve the tissue-specific gene networks:
* [HumanBase](https://hb.flatironinstitute.org)

To retrieve synthetic lethality information: 
* [SynLethDB 2.0](https://synlethdb.sist.shanghaitech.edu.cn/v2/)

The training data of this project can be accessed from the public repository Open Science Framework (OSF):
* DDR combination treatment responses and dose-response matrix: https://osf.io/8hbsx/
* molecular signatures: https://osf.io/8mxgj/


## Repository Structure

This repository contains the following code, where all the rest of data in this project can be carried out in the following steps

### Data Preparation

- `./feature`: before we run the model, all required data, QC and features are preprocessed in this folder. it includes the following subdirectories:

  - `./data`
  - `./QC`
  - `./chemical_structure_features`
  - `./drug_similarity_features`  
  - `./mode_of_action`
  - `./monotherapy_features`
  - `./target_gene`
  - `./geneset_features`
  - `./molecular_features`
  - `./synthetic_lethality_features`
  - `./tissue_specific_networks`

From the above subdirectories:

- `data`: deposit the drug response data, geneset, chemical structure, molecular readouts etc.
 
  it contains the main training dataset:
  - `merck_confidential_responses_training_20200221.tsv`: training and testing sets, including combination therapy and monotherapy. the later will be used as features.
    
  - `merck_confidential_genesets_20200414.tsv`: geneset annotations from literature
    
  - `merck_confidential_features_20200221.tsv`: molecular features for all cell lines on gene level, including exp, lof, snv, cnv, lof_pat, coh_pat, ddr

  the hold-out validation dataset is also under this folder:
  - `./hold_out`: hold out validation dataset
      
    - `./ho1`: Hold-out set 1 asks for responses of previously unobserved cell lines from known tissue types. 
    - `./ho2`: Hold-out set 2, you have not seen any responses for these tissues yet, so your model will be confronted with new entries in the ".meta_cancer_*" features that it hasn't seen before. 

All the following subdirectories in the `./features` folder are for storage of different features of the model:

- `chemical_structure_features`
- `drug_similarity_features`  
- `mode_of_action`
- `monotherapy_features`
- `target_gene`
- `geneset_features`
- `molecular_features`
- `synthetic_lethality_features`
- `tissue_specific_networks`

To preprocess and quality check of all information above, we use code from the folder `QC`. it contains the following programs:

- `QC.py`: generate all features
  
  run

  ```
  python QC.py
  ```
    
  it contains the following modules that can generate following features:
    
  - geneset_qc: 
    - `./geneset_features/`

  - molecular_qc: 
    
    - `./molecular_features/`
    - `./QC/ddr_genes.txt`
    
  - response_qc: 
    - `./monotherapy_features/monotherapy_features.tsv`
    - `./reproducibility/`
    - `./target_genes/all_drugs_summary.csv`

- `QC_visualization.ipynb`: visualize the reproducibility, number and quality of all datasets
- `pull_drug_target.py`: pull drug target information from public databases we mentioned in the [Dependencies](./README.md##Dependencies).
  
  it will execute following tasks:
    - `./target_genes/`
      - update `all_drugs_summary.csv`
      - update `drug2gene.json`
    
- `split_train_test.py`: split the training and testing dataset for the k-fold cross validation, using files in the `./data` folder we mentioned above. the training and testing dataset can be slipted by cell line/tissue types/drug/moa etc. Running this program will generate the following folders in the parent directory:
    - `test_by_cell_line`
      - `fold_0`
        - `Train.tsv`
        - `Test.tsv`
      - ... ...
      - `fold_9`
        - `Train.tsv`
        - `Test.tsv`
      - `hold_out_validation`
        - `Test.tsv`


    - `test_by_indication  
      - `fold_bladder`
        - `Train.tsv`
        - `Test.tsv`
      - ... ...
      - `fold_sarcoma`
        - `Train.tsv`
        - `Test.tsv`
      - `hold_out_validation`
        - `Test_cervix.tsv`
        - `Test_kidney.tsv`


### Train Final Models and Create SHAP Analysis Plots

- `./master`: the code for model training and evaludation
  - `./master_code`: LightGBM model for response prediction, which consists of the following python programs and subdirectory:
    - `./main.py`
    - `./build_feature_dataset.py`
    - `./common.py`
    - `./models.py`
    - `./shap_analysis.py`
    - `./utils.py`
    - `./downstrean_analysis/`

To check the usage of the program, we run the following command inside the `./master_code` direcotry:

```
python main.py --help
```
it would show the instructions on running this model:

```
usage: main.py [-h] [-p PATH] [--exclude_synergy_batch] [--exclude_cell_line] [--exclude_cancer_type] [-s PREDICTED_SCORE] [-f FEATURES [FEATURES ...]]
               [--mol_type MOL_TYPE] [--surrogate_gene SURROGATE_GENE] [--surrogate_geneset SURROGATE_GENESET] [--surrogate_synleth SURROGATE_SYNLETH]
               [--surrogate_chem SURROGATE_CHEM] [--hold_out_test] [--evaluate_shap EVALUATE_SHAP [EVALUATE_SHAP ...]] [--skip_training]

Build drug synergy prediction machine learning models.

options:
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
                                    network: tissue-specific gene network information from HumanBase.
                                    geneset: 9,450 geneset features;
                                    chemical_structure: drug chemical structure fingerprints;
                                    synthetic_lethality: synthetic lethality information
                                    (default: [monotherapy, moa, drug_name, molecular, target_gene, geneset, chemical_structure, synthetic_lethality])
                                    
  --mol_type MOL_TYPE    number of top genes for surrogate model:
                                    chose one from: exp, cnv, snv, lof, coh_pat, lof_pat, ddr
                                    if not specified, default as None
                                    
  --surrogate_gene SURROGATE_GENE
                         number of top genes for surrogate model:
                                    chose one from: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 250, 500, 1000, 2000
                                    if not specified, default as None
                                    
  --surrogate_geneset SURROGATE_GENESET
                         number of top genes for surrogate model:
                                    chose one from: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200
                                    if not specified, default as None
                                    
  --surrogate_synleth SURROGATE_SYNLETH
                         number of top synthetic lethal pairs for surrogate model:
                                    chose one from: 5, 10, 15, 20, 25, 30, 40, 50, 80, 100, 200, 500
                                    if not specified, default as None
                                    
  --surrogate_chem SURROGATE_CHEM
                         type of chemical stucture feature fingerprint for surroagte model:
                                    chose one from: MACCS, Morgan, FP2, FP3, FP4, RDK
                                    if not specified, default as None
                                    
  --hold_out_test       if used, test on the hold-out set
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

There are sample code to run the master code program using different criterias in the parent directory:

- `../bash.sh`
- `../bash_holdout.sh`
- `../bash_mol.sh`
- `../bash_shap.sh`
- `../bash_topgenes.sh`
- `../bash_topsynleth.sh`

### Downstream SHAP feature importance analysis

This part is crucial for analysing potential molecualr biomarkers. 

First, in the following directory:
```
cd ./downstream_analysis
```

run the following code:
```
python statistical_analysis.py
``` 

This program will carry out analysis in the following aspects:

* 1. Analyse SHAP results by feature category.
* 2. Analyse ddr gene's feature importance by it's molecular features (cnv, snv, lof and exp)
* 3. Analyse tissue-specific ddr gene features.
    

### Result Summarization

Check the following directory for the result summarization"

```
cd ./result_summary
```
the following programs are in this folder:

- `calculate_efficacy.ipynb`: generate the treatment prioritizatiion results, compared to the basedline treatment. (in different cell lines/treatments/moas)
- `generate_surrogate_model.ipynb`: generate surrogate model features based on shap analysis of the full model
- `result_summary.py`: summarize the prediction results in cross-validation and hold-out validation in all different validation settings
- `tissue-specific_results.py`: the prediction results on different tissue types (of the full model)
- `utils.py`: utilities for the above code

### Post Analysis

Other statistical and visualization results in this paper.

```
cd ./molceular_analysis
```
the following programs are in this folder:
- `gene_interaction_synergy.py`: generate top interacted genes with the core DDR targets (ATR/ATM/DNAPK) for efficacy and synergy. 
- ` monotherpay_top_genes.py`: top genes positively/negatively correlated with efficacy in ATRi monotherapy.ge  


## Reference

TBD

## Manuscript preparation
The manuscript/figure/data in this paper are organized by following the format requirements of submissions to Nature Communications.
- `./Code for figures`
  - `./plot.R`: generate figures using data
  - `./data`: pulled directly from project folder
  - `./souce_data`: the source data for all figures (required by Nature Communications)
  - `./Figures`: all main and supp figures



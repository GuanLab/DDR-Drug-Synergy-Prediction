#!usr/bin/python3
#author:@rayezh
import pandas as pd
import numpy as np
import argparse, textwrap, sys
from common import *

def main():
    
    parser = argparse.ArgumentParser(description = "Build drug synergy prediction machine learning models.",
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-p','--path', type = str,
            help = "Path to cross-validation split data",
            default = '../../test_by_cell_line')
    
    parser.add_argument('--exclude_synergy_batch', 
            action = 'store_true',
            help = '''If specified, exclude synergy batch as feature. If path is 'test_by_batch', exclude synergy batch automatically.''')
    
    parser.add_argument('--exclude_cell_line', 
            action = 'store_true',
            help = '''If specified, exclude cell line as feature. If path is 'test_by_cell_line', exclude cell line automatically.''')
    
    parser.add_argument('--exclude_cancer_type', 
            action = 'store_true',
            help = '''If specified, exclude cancer type/subtype as feature. If path is 'test_by_indication', exclude cancer type automatically.''')
    
    parser.add_argument('-s', '--predicted_score', type = str,
            help = '''Synergy measurement as training target.
            aoc: drug efficacy;
            bliss: drug synergy;
            (default: aoc)
            ''',
            default = 'aoc')
    
    parser.add_argument('-f', '--features', nargs = '+',
            type = str,
            default =['monotherapy', 'moa', 'drug_name','molecular', 'target_gene', 'geneset', 'chemical_structure'],
            help = '''Feature options in building feature set:
            monotherapy: monotherapy responses;
            moa: mode-of-action categorical feature;
            drug_name: drug name as categorical feature;
            molecular: four types of single gene molecular biomarkers,two types of gene cluster biomarkers,78 ddr markers;
            target_gene: target gene information (target gene information can only be used when molecular information is included);
            geneset: 9,450 geneset features;
            chemical_structure: drug chemical structure fingerprints;
            (default: [monotherapy, moa, drug_name, molecular, target_gene, geneset, chemical_structure, dependency, synthetic_lethality])
            ''')
    
    parser.add_argument('--surrogate_gene', type = str,
            required = False,
            help = ''' number of top genes for surrogate model:
            chose one from: 2000, 1000, 500, 250, 200, 125, 100, 75, 50, 25, 20, 10,
            if not specified, default as None
            '''
            )
    
    parser.add_argument('--surrogate_geneset', type = str,
            required = False,
            help = ''' number of top genes for surrogate model:
            chose one from: 2000, 1000, 500, 250, 200, 125, 100, 75, 50, 25, 20, 10,
            if not specified, default as None
            '''
            )

    parser.add_argument('--surrogate_chem', type = str,
            required = False,
            help = ''' type of chemical stucture feature fingerprint for surroagte model:
            chose one from: MACCS, Morgan, FP2, FP3, FP4, RDK
            if not specified, default as None
            ''')
    parser.add_argument('--held_out_test', action = 'store_true',
            help = 'if used, test on the held-out set')
    parser.add_argument('--evaluate_shap', nargs = '+',
            type = str,
            required=False,
            help = '''subset type for subset specific SHAP analysis.
            example:
            'Train'
            'Test'
            'moa': Test split by ATM/ATR/DNAPK
            'all_moa': Test split by moa pairs
            'cell_line'
            'tissue'
            ''')
    parser.add_argument('--skip_training', action = 'store_true',
            help = 'if used, skip training.')
    args = parser.parse_args()
    print(args.features)
    
    # if test by batch: exclude batch as synergy feature
    if 'test_by_batch' in args.path:
        print("Cross-validation split by batch: exclude synergy batch feature!")
        args.exclude_synergy_batch = True
    if 'test_by_cell_line' in args.path:
        print("Cross-validation split by cell line: exclude cell line feature!")
        args.exclude_cell_line = True
    if 'test_by_indication' in args.path:
        print("Cross-validation split by indication: exclude cell line and cancer type feature!")
        args.exclude_cancer_type = True
        args.exclude_cell_line = True
    if ('molecular' not in args.features) and('target_gene' in args.features):
        sys.exit("Can not use target gene information without molecular feature!")
    if args.surrogate_gene is not None:
        if args.surrogate_gene not in ['2000', '1000', '500', '250', '200', '125', '100', '75', '50', '25', '20', '10']:
            sys.exit('Surrogate gene model not available!')
        print('Train surrogate model with top '+args.surrogate_gene+' genes based on SHAP analysis.')
    if args.surrogate_geneset is not None:
        if args.surrogate_geneset not in ['2000', '1000', '500', '250', '200', '125', '100', '75', '50', '25', '20', '10']:
            sys.exit('Surrogate geneset model not available!')
        print('Train surrogate model with top '+args.surrogate_geneset+' genesets based on SHAP analysis.')
    if args.surrogate_chem is not None:
        if args.surrogate_chem not in ['MACCS', 'Morgan', 'FP2', 'FP3', 'FP4', 'RDK']:
             sys.exit('Surrogate chem model not available!')
        print('Train surrogate model with '+args.surrogate_chem+' chemical structure fingerprints.')
    if args.held_out_test == True:
        print("Will test on held-out set")
    if args.evaluate_shap is not None:
        print(args.evaluate_shap)
    if args.skip_training == True:
        print("Skip training ...")
    opts = vars(args)
    run(**opts)


def run(path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, surrogate_gene, surrogate_geneset, surrogate_chem, held_out_test, evaluate_shap, skip_training):
    """
    """
    if skip_training:
        pass
    else:
        eva_df, eva_all = five_fold_cross_validation(path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, surrogate_gene, surrogate_geneset, surrogate_chem)
        eva_df["features"] = "+".join(features)
        eva_df.to_csv("cor_"+predicted_score+".csv", index = False)
        result = open("cv-95ci_"+predicted_score+".txt", 'w')
        result.write(eva_all)
        result.close()
    
    #if held_out_test:
    #    if path == '../../test_by_cell_line':
    #        hold_out_path = '../../held_out_stage_2/ho1/synergy_responses_with_monotherapy.tsv'
    #        predict_on_held_out(hold_out_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, surrogate_gene, surrogate_geneset, surrogate_chem)
    #    if path == '../../test_by_indication':
    #        hold_out_path = '../../held_out_stage_2/ho2/synergy_responses_with_monotherapy.tsv'
    #        predict_on_held_out(hold_out_path,exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, surrogate_gene, surrogate_geneset, surrogate_chem)

    # carry out SHAP
    if evaluate_shap is not None:
        for subset_type in evaluate_shap:
            evaluate_SHAP_analysis(path, subset_type, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, surrogate_gene, surrogate_geneset, surrogate_chem)

if __name__ == "__main__":
    main()




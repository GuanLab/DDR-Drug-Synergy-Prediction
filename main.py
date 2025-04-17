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
            default =['monotherapy', 'moa', 'drug_name','molecular', 'target_gene', 'network', 'geneset', 'chemical_structure', 'synthetic_lethality'],
            help = '''Feature options in building feature set:
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
            ''')
    
    parser.add_argument('--mol_type', type = str,
            required = False,
            help = ''' number of top genes for surrogate model:
            chose one from: exp, cnv, snv, lof, coh_pat, lof_pat, ddr
            if not specified, default as None
            '''
            )

    
    parser.add_argument('--surrogate_gene', type = str,
            required = False,
            help = ''' number of top genes for surrogate model:
            chose one from: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 250, 500, 1000, 2000
            if not specified, default as None
            '''
            )
    
    parser.add_argument('--surrogate_geneset', type = str,
            required = False,
            help = ''' number of top genes for surrogate model:
            chose one from: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200
            if not specified, default as None
            '''
            )
    
    parser.add_argument('--surrogate_synleth', type = str,
            required = False,
            help = ''' number of top synthetic lethal pairs for surrogate model:
            chose one from: 5, 10, 15, 20, 25, 30, 40, 50, 80, 100, 200, 500
            if not specified, default as None
            ''')

    parser.add_argument('--surrogate_chem', type = str,
            required = False,
            help = ''' type of chemical stucture feature fingerprint for surroagte model:
            chose one from: MACCS, Morgan, FP2, FP3, FP4, RDK
            if not specified, default as None
            ''')
    
    parser.add_argument('--hold_out_test', action = 'store_true',
            help = 'if used, test on the hold-out set')
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
    if ('molecular' not in args.features)and('synthetic_lethality' not in args.features)and('target_gene' in args.features):
        sys.exit("Can not use target gene information without molecular feature!")
    if ('molecular' not in args.features) and('mol_type' in args.features):
        sys.exit("Can not use molecular subtype information without molecular feature!")
    if args.surrogate_gene is not None:
        if args.surrogate_gene not in ['10', '20', '30', '40', '50', '60', '70', '80', '90', '100', '150', '250', '500', '1000', '2000']:
            sys.exit('Surrogate gene model not available!')
        print('Train surrogate model with top '+args.surrogate_gene+' genes based on SHAP analysis.')
    if args.surrogate_geneset is not None:
        if args.surrogate_geneset not in ['2000', '1000', '500', '250', '200', '125', '100', '75', '50', '25', '20', '10']:
            sys.exit('Surrogate geneset model not available!')
        print('Train surrogate model with top '+args.surrogate_geneset+' genesets based on SHAP analysis.')
    if args.surrogate_synleth is not None:
        if args.surrogate_synleth not in ['5', '10', '15', '20', '25', '30', '40', '50', '80', '100', '200', '500']:
            sys.exit('Surrogate synthetic lethality model not available!')
        print('Train surrogate model with top '+args.surrogate_synleth+' synthetic lethal pairs based on SHAP analysis.')
    if args.surrogate_chem is not None:
        if args.surrogate_chem not in ['MACCS', 'Morgan', 'FP2', 'FP3', 'FP4', 'RDK']:
             sys.exit('Surrogate chem model not available!')
        print('Train surrogate model with '+args.surrogate_chem+' chemical structure fingerprints.')
    if args.hold_out_test == True:
        print("Will test on hold-out set")
    if args.evaluate_shap is not None:
        print(args.evaluate_shap)
    if args.skip_training == True:
        print("Skip training ...")
    opts = vars(args)
    run(**opts)


def run(path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, mol_type, surrogate_gene, surrogate_geneset, surrogate_synleth, surrogate_chem, hold_out_test, evaluate_shap, skip_training):
    """
    """
    from glob import glob
    if skip_training:
        pass
    else:
        eva_df, eva_all = k_fold_cross_validation(path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, mol_type, surrogate_gene, surrogate_geneset, surrogate_synleth, surrogate_chem)
        eva_df["features"] = "+".join(features)
        eva_df.to_csv("cor_"+predicted_score+".csv", index = False)
        result = open("cv-95ci_"+predicted_score+".txt", 'w')
        result.write(eva_all)
        result.close()
    
    if hold_out_test: # test the best model on hold-out set
        hold_out_path = path+'/hold_out_validation/*.tsv'
        os.makedirs('./hold_out_validation', exist_ok = True)
        for f in glob(hold_out_path):
            subtype = f.split('/')[-1].split('.')[0]
            print(subtype) # for cell line, print Test; for indication, print Test_kidney and Test_cervix
            new_path = './hold_out_validation/'+subtype+'/'
            os.makedirs(new_path, exist_ok = True)

            # evaluate SHAP analysis
            evaluate_SHAP_analysis_on_holdout(path, 'Test', exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, mol_type, surrogate_gene, surrogate_geneset, surrogate_synleth, surrogate_chem)

            eva_df, eva_all = predict_on_hold_out(f, new_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, mol_type, surrogate_gene, surrogate_geneset, surrogate_synleth, surrogate_chem)
            eva_df["features"] = "+".join(features)
            eva_df.to_csv(new_path+'cor_'+predicted_score+'.csv', index = False)
            result = open(new_path+'cv-95ci_'+predicted_score+'.txt', 'w')
            result.write(eva_all)
            result.close()
            
    # carry out SHAP
    if evaluate_shap is not None:
        for subset_type in evaluate_shap:
            evaluate_SHAP_analysis(path, subset_type, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, predicted_score, features, mol_type, surrogate_gene, surrogate_geneset, surrogate_synleth, surrogate_chem)

if __name__ == "__main__":
    main()




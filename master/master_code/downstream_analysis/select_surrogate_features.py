import pandas as pd
import pickle
import os

"""
Select features for surrogate model 

prerequisite:
    ../features_molecular_target_gene_network/aoc[bliss]_*.tsv

"""

def generate_surrogate_list(m_version):
    """
    params
    ------
    m_version: model version (geneset or gene)
    """
    # model version to select feature from
    if m_version == 'gene':
        model_path = '../../features_molecular_target_gene_network/'
        feature_type = "Molecular_exp" # feature type to select from
    if m_version == 'geneset':
        model_path = '../../features_geneset/'
        feature_type = "Geneset"

    for score in ['aoc', 'bliss']:
        df = pd.read_csv(model_path+score+'_totalSHAP.tsv', sep = '\t')
        df_sub = df[df['feature_type'] == feature_type].copy()

        # select new top n features (currently just molecular)
        for n in [2000, 1000, 750, 500, 250, 200, 125, 100, 75, 50, 25, 20, 10]:
            new_features = df_sub.groupby('feature', as_index=False).agg('mean').sort_values(by='SHAP_val', ascending = False).reset_index(drop = True)['feature'].iloc[:n]
            #print(new_features)
            pickle.dump(new_features, open('./surrogate_models/'+score+'_surrogate_'+m_version+'_'+str(n)+'.pkl', 'wb'))

os.makedirs('./surrogate_models', exist_ok = True)
generate_surrogate_list('gene')
generate_surrogate_list('geneset')


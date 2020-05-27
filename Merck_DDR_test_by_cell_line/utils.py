#!usr/bin/python3
#author:@rayezh

def pd_to_XY(target, data, if_train = False):
    import pandas as pd
    import numpy as np
    import json, glob,os
    from tqdm import tqdm

    target_genes = []
    #read moa's target genes
    f = open('../data/moa_to_genes.json','r').read()
    moa2gene = json.loads(f)
    
    gene_features = pd.read_csv('../data/merck_confidential_features_20200221.tsv', sep ='\t', index_col = 0)
    cols = gene_features.columns.to_list()
    embed = {}
    for i in range(13):
        embed.update({'Oncolead_'+'{:0>3}'.format(i):np.array([i])})
    df = []
    mono_cols = open('feature_header.txt','r').read().rstrip().split(',')
    feature_cols = ['Synergy_batch']+mono_cols+gene_features.columns[3:].to_list()+['.response_'+target]
    print(len(feature_cols))
    for _, row in tqdm(data.iterrows(), total = data.shape[0]):
        try:
            feature0 =embed[row['.identifier_batch']]
            feature1 = np.loadtxt('./cell_drug_features/'+row['.identifier_sample_name']+'.'+row['.metadata_treatment_1'])
            feature2 = np.loadtxt('./cell_drug_features/'+row['.identifier_sample_name']+'.'+row['.metadata_treatment_2'])
            target_genes = list(set(moa2gene[row['.metadata_moa_1']]+moa2gene[row['.metadata_moa_2']])) #overlapped target genes of two drugs
            cols = gene_features.columns.to_list()
            tmp = [subs+'_exp' for subs in target_genes]
            tmp = list(set(tmp)&set(cols))
            target_genes = tmp
            mol = gene_features.loc[gene_features['.identifier_sample_name'] == row['.identifier_sample_name']].copy()
            mol[target_genes] = 0
            feature3 = np.array(mol.iloc[0,3:].to_list())
            #print(feature3.shape)
            y = np.array([row['.response_'+target]])
            x = np.concatenate((feature0,feature1,feature2, feature3, y))
            #print(x.shape)
            df.append(x)
            if if_train:
                x = np.concatenate((feature0,feature2,feature1, feature3, y))
                df.append(x)
            #Y.append(y)
        except:
            pass
    df= np.array(df)
    #print(len(feature_cols))
    #print(X.shape)
    df = pd.DataFrame(df, columns=feature_cols)
    X = df.iloc[:,:-1]
    Y = df['.response_'+target]
    return X, Y


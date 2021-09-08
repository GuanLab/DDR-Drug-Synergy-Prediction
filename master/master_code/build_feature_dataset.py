#usr/bin/env python3
#author: @rayezh

import pandas as pd
import numpy as np
import json, glob, os, re, pickle
from collections import defaultdict
from tqdm import tqdm
from utils import encode_categorical_feature

# use average monotherapy instead of separated features

def build_feature_dataset(data, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem, if_train = True, return_label = True):
    """ Construc feature set
    
    Parameters:
    -----------
    
    data: Pandas dataframe
        synergy response data in training or testing dataset
    
    exclude_synergy_batch: boolean
        if True, exclude synergy batch.
        Used for cross-batch validation

    exclude_cell_line: boolean
        if True, exclude cell line
        Used for cross-cell line validation

    exclude_cancer: boolean
        if True, exclude cancer tissue
        Used for cross-indication validation

    target: string
        prediction target('aoc' or 'bliss')
    
    features: list
        features to use in prediction.
        monotherapy response
        molecluar
        chemical_structure
        crispr

    if_train: boolean
        if construct training data, true; if test data, false
    
    return_label: boolean
        option to return label (for cross validation)
        or not to return label (for held-out set prediction, where the label is unknown))
    
    Yields:
    -------
    feature_name: list of strings
        feature names
    X: Numpy array
        feature matrix
    Y: Numpy array
        predicted label(aoc score or bliss score)
    """

    # load feature datasets
    # monotherapy
    if 'monotherapy' in features:
        all_monotherapy = pd.read_csv('../../feature/monotherapy_features/monotherapy_features.tsv', sep = '\t', index_col = 0)
    # chemical structure
    if 'chemical_structure' in features:
        all_chemical_structure = pd.read_csv('../../feature/chemical_structure_features/all_chemical_structure.csv', index_col = 0)
        if surrogate_chem is not None:
            cols = [i for i in all_chemical_structure.columns if i.startswith(surrogate_chem)]
            all_chemical_structure = all_chemical_structure[cols]
            #print(all_chemical_structure)

    # molecular markers
    if 'molecular' in features:
        all_molecular = pd.read_csv('../../feature/molecular_features/quantiled_molecular_features.tsv', sep = '\t', index_col = 0)
        # TODO: add mol_type = ['exp', 'cnv', 'snv', 'lof', 'coh_pat', 'lof_pat', 'ddr'] 
        #all_molecular = pd.read_csv('../../feature/molecular_features/'+mol_type+'_molecular_features.tsv', sep = '\t', index_col = 0)  # coh_pat
        if surrogate_gene is not None:
            surrogate_features = pickle.load(open('./downstream_analysis/surrogate_models/'+target+'_surrogate_gene_'+str(surrogate_gene)+'.pkl', 'rb'))
            all_molecular = all_molecular[surrogate_features]
            #print(all_molecular)
    
    # geneset annotations
    if 'geneset' in features:
        geneset = pd.read_csv('../../feature/geneset_features/geneset_features.csv', index_col = 0)
        if surrogate_geneset is not None:
            surrogate_features = pickle.load(open('./downstream_analysis/surrogate_models/'+target+'_surrogate_geneset_'+str(surrogate_geneset)+'.pkl', 'rb'))
            surrogate_features =  [x.split('Geneset_')[1] for x in surrogate_features]
            geneset = geneset.loc[list(surrogate_features),]
            #print(geneset)

    if 'dependency' in features:
        # depmap dependency
        ceres_cancer_dependency = json.load(open('../../feature/cancer_dependency_features/integrated_Sanger_Broad_essentiality_matrices_20200402/CERES_FC_dep.json', 'r'))
        crisprclear_cancer_dependency = json.load(open('../../feature/cancer_dependency_features/integrated_Sanger_Broad_essentiality_matrices_20200402/CRISPRcleanR_FC_dep.json', 'r'))

    # synthetic lethality
    if 'synthetic_lethality' in features:
        synleth = pd.read_csv('../../feature/synthetic_lethality_features/new_Human_SL.csv')
        synleth_dic = {}
        for _,r in synleth.iterrows(): 
            synleth_dic.update({(r['gene_a.name'],r['gene_b.name']):r["SL.statistic_score"]})
            #synleth_dic.update({(r['gene_b.name'],r['gene_a.name']):r["SL.statistic_score"]})
        #gene rank
        exp2rk = pickle.load(open('../../feature/QC/exp_rank.pkl', 'rb'))

    if 'target_gene' in features:
        # load drug-specific target gene information
        drug2gene = json.load(open('../../feature/target_gene/drug2gene.json','r'))
        drug2gene = defaultdict(lambda:[], drug2gene)

    # load gene network information from mousenet
    if 'network' in features:
        network = pickle.load(open('../../feature/mousenet_features/target_gene_network.pkl','rb'))

    # categotical encoding
    if 'moa' in features:
        encode_moa = json.load(open('../../feature/categorical_embedding_features/moa.json', 'r'))
    if 'drug_name' in features:
        encode_treatment = json.load(open('../../feature/categorical_embedding_features/treatment.json', 'r'))
    if not exclude_synergy_batch:
        encode_batch = json.load(open('../../feature/categorical_embedding_features/batch.json', 'r'))
    if not exclude_cell_line:
       encode_cell_line = json.load(open('../../feature/categorical_embedding_features/cell_line.json', 'r')) 
    if not exclude_cancer_type:
        encode_cancer_type = json.load(open('../../feature/categorical_embedding_features/cancer_type.json', 'r'))
        encode_cancer_subtype = json.load(open('../../feature/categorical_embedding_features/cancer_subtype.json', 'r'))

    def build_features(row, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, features):
        """ Build up all types of features
        
        Params
        ------
        row: Pandas row series
            a single line of synergy experimental data
        features: a list of str
            features used in feature set
        exclude_synergy_batch: boolean
        exclude_cell_line: boolean
        exclude_cancer_type: boolean

        Yields:
        -------
        feature_names: list
        feature_values: list
        """

        def synergy_feature(row, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type):
            """ Synergy feature
            .identifier_batch       
            .identifier_sample_name
            .metadata_cancer_type   
            .metadata_cancer_subtype
            .metadata_treatment_remarks

            Params
            ------
            row: Pandas row series
            exclude_synergy_batch: boolean
            exclude_cell_line: boolean
            exclude_cancer_type: boolean
            
            Yields
            ------
            feature_names: list
            feature_values: list
            """
            features = {}

            if exclude_synergy_batch:
                pass
            else:
                try:
                    features.update({"Synergy_batch": encode_batch[row['.identifier_batch']]})
                except:
                    features.update({"Synergy_batch": np.nan})
            
            if exclude_cell_line:
                pass
            else:
                try:
                    features.update({"Cell_line": encode_cell_line[row['.identifier_sample_name']]})
                except:
                    features.update({"Cell_line": np.nan})
            
            if exclude_cancer_type:
                pass
            else:
                try:
                    features.update({"Cancer_type": encode_cancer_type[row['.metadata_cancer_type']]})
                except:
                    features.update({"Cancer_type": np.nan})
                try:
                    features.update({"Cancer_subtype": encode_cancer_subtype[row['.metadata_cancer_subtype']]})
                except:
                    features.update({"Cancer_subtype": np.nan})
            """
            try:
                features.update({"Concentration": encode_remark[str(row['.metadata_treatment_remarks'])]})
            except:
                features.update({"Concentration": np.nan})
            """
            feature_names = list(features.keys())
            feature_values = list(features.values())
            
            return feature_names, feature_values
        
        def moa_feature(row):
            """ Mode-of-action as categorical feature

            Params
            ------
            row: Pandas row series

            Yields
            ------
            feature_names: list
            feature_values: list
            """
            feature_names = []
            feature_values = []
            for i in ['1','2']:
                moa = row['.metadata_moa_'+i]
                feature_names.extend(['Treatment_'+i+'_moa'])
                feature_values.extend([encode_moa[moa]])

            return feature_names, feature_values

        def drug_name_feature(row):
            """ Drug name as categorical feature
            
            Params
            ------
            row:Pandas row series

            Yields
            ------
            feature_names: list
            feature_values: list
            """
            feature_names = []
            feature_values = []
            for i in ['1','2']:
                treatment = row['.metadata_treatment_'+i]
                feature_names.extend(['Treatment_'+i+'_name'])
                feature_values.extend([encode_treatment[treatment]])

            return feature_names, feature_values

        def monotherapy_feature(row):
            """ Monotherapy response feature
            
            Params
            ------
            row: Pandas row series
            
            Yields
            ------
            feature_names: list
            feature_values: list
            """
            feature_names = []
            feature_values = []
            
            for i in ['1','2']:
                cell_line_treatment = row['.identifier_sample_name']+'.'+row['.metadata_treatment_'+i]
                #feature_names.extend(['Treatment_'+i+'_'+b for b in all_monotherapy.columns])
                feature_names.extend(['Treatment_'+i+'_ave'])    # one average monotherapy response from all batches
                if cell_line_treatment in all_monotherapy.index:
                    #feature_values.extend(all_monotherapy.loc[cell_line_treatment].to_list())
                    feature_values.extend([np.nanmean(all_monotherapy.loc[cell_line_treatment])]) # average response from all batches 
                else:
                    #feature_values.extend([np.nan]*len(all_monotherapy.columns))
                    feature_values.extend([np.nan]) # replace average with nan
                    print("No monotherapy response data of "+ cell_line_treatment+"! Replace with NaN!")
                    
            return feature_names, feature_values

        def molecular_feature(row, include_target_gene = False, include_network =  False):
            """ Molecular features
            
            Params
            ------
            row: Pandas row series
            include_target_gene: boolean
                if include target gene information or not
            include_network: boolean
                include network information or not
                can only be set to true if include_target_gene is true

            Yields
            ------
            feature_names: list
            feature_values: list
            """

            mol = all_molecular.loc[row['.identifier_sample_name']].copy()
            
            # include target gene information
            if include_target_gene:
                # get target genes from treatment 1 and 2
                target_genes = set(drug2gene[row['.metadata_treatment_1']]+drug2gene[row['.metadata_treatment_2']])
                target_genes_exp = [g+'_exp' for g in target_genes]
                # set <target gene>_exp to 0
                for exp in target_genes_exp:
                    if exp in mol.index:
                        mol[exp] = 0
                    else:
                        pass

                if include_network:
                    # maximum relationship with target genes
                    p2target = defaultdict(lambda:0)
                    for target_gene in target_genes:
                        for k, v in network[target_gene].items():
                            if v>p2target[k]:
                                p2target[k] = v
                    for k,v in p2target.items():
                        if k+'_exp' in mol.index:
                            mol[k+'_exp'] = mol[k+'_exp']*(1-v)

            feature_names = list(mol.index)
            feature_values = list(mol)
            
            return feature_names, feature_values

        def chemical_structure_feature(row):
            """ Drug features:
                1. treatment name (categorical feature)
                2. drug chemical structure feature
            
            Params
            ------
            row: Pandas row series

            Yields
            ------
            feature_names: list
            feature_values: list
            """
            feature_names = []
            feature_values = []
            for i in ['1','2']:
                treatment = row['.metadata_treatment_'+i]
                # extend feature name
                feature_names.extend(['Treatment_'+i+'_'+c for c in all_chemical_structure.columns.to_list()])
                # extend feature values
                if treatment in all_chemical_structure.index:
                    feature_values.extend(all_chemical_structure.loc[treatment].to_list())
                else:
                    feature_values.extend([np.nan]*len(all_chemical_structure.columns))
                    if treatment == "GAMMA":
                        pass
                    else:
                        print("No chemical structure data of "+treatment+"! Replace with NaN!")
            
            return feature_names, feature_values
        
        def geneset_feature(row):
            """ Geneset feature 
            
            Params
            ------

            Yields
            ------
            feature_value: a vector of n genesets by 1.
                for each geneset, it returns a number denoting how many genes targeted by treatments are in the pathway
            """
            feature_names = list(geneset.index)
            feature_values = np.zeros(len(feature_names))

            # get target genes from treatment 1 and 2
            target_genes = set(drug2gene[row['.metadata_treatment_1']]+drug2gene[row['.metadata_treatment_2']])

            feature_names = ["Geneset_"+x for x in list(geneset.index)]
            for g in target_genes:
                try:
                    feature_values += np.array(geneset[g])
                except:
                    pass
                    #print("no geneset feature!")
            feature_values = list(feature_values)
            
            return feature_names, feature_values

        def cancer_dependency_feature(row):
            """ Cancer dependency feature
            
            Params
            ------
            row: Pandas row series
            
            Yields
            ------
            feature_names: list
            feature_values: list
            """
            cell = row['.identifier_sample_name']
            moas = [row['.metadata_moa_1'], row['.metadata_moa_2']]
            
            feature_names = []
            feature_values = []

            for i in ['ceres', 'crispr']:
                f_name = 'CRISPR_cancer_dependency_'+i
                if i == 'ceres':
                    df = ceres_cancer_dependency
                if i == 'crispr':
                    df = crisprclear_cancer_dependency
                for j in ['max']: # max dependency
                    f_name = f_name+'_'+j
                    feature_names.append(f_name)
                    try:
                        v = []
                        for moa in moas:
                            v.append(df[moa][cell][j])
                        if j == 'max':
                            v = max(v)
                        elif j == 'mean':
                            v = np.mean(v)
                        elif j == 'min':
                            v = min(v)
                    except:
                        v = np.nan
                    feature_values.append(v)

            return feature_names, feature_values

        def synthetic_lethality_feature(mol_name, mol_feature):
            """ Synthetic lethality feature

            Params
            ------
            row: Pandas row series

            Yields
            ------
            feature_names: list
            feature_values: list

            """
            feature_names = []
            feature_values = []

            #exp
            idx_exp = [i for i,j in enumerate(mol_name) if j.split('_')[1] =='exp']
            f_exp = {mol_name[i].split('_')[0]:mol_feature[i] for i in idx_exp}
            n_exp = [mol_name[i].split('_')[0] for i in idx_exp]

            #lof
            idx_lof = [i for i,j in enumerate(mol_name) if j.split('_')[1] =='lof']
            f_lof = {mol_name[i].split('_')[0]:mol_feature[i] for i in idx_lof}
            n_lof = [mol_name[i].split('_')[0] for i in idx_lof]

            feature_exp = 0
            #feature_lof = 0
            
            for k in synleth_dic:
                gene_a = k[0]
                gene_b = k[1]
                syn = synleth_dic[k]
                if (gene_a in n_exp) and (gene_b in n_exp):
                    feature_n = "Synleth_"+('_').join(sorted([gene_a,gene_b]))
                    #print(gene_a, np.round(f_exp[gene_a],4))
                    #print(gene_b, np.round(f_exp[gene_b],4))
                    #feature_v = max(exp2rk[gene_a][np.round(f_exp[gene_a],4)],exp2rk[gene_b][np.round(f_exp[gene_b],4)])
                    rk_a = f_exp[gene_a]/exp2rk[gene_a]['median']
                    rk_b = f_exp[gene_b]/exp2rk[gene_b]['median']
                    #feature_v = max(rk_a,rk_b)
                    feature_v = 1/(rk_a*rk_b+0.001)
                    feature_names.append(feature_n)
                    feature_values.append(feature_v)
                #if (gene_a in n_lof) and (gene_b in n_lof):
                #    lof = 1-syn*f_lof[gene_a]*f_lof[gene_b]
                #    if feature_lof < lof:
                #        feature_lof = lof

            return feature_names, feature_values

        # aggregate all features above:
        feature_names = []
        feature_values = []
        
        syn_name, syn_feature = synergy_feature(row, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type)
        feature_names += syn_name
        feature_values += syn_feature

        if 'moa' in features:
            moa_name, moa_feature = moa_feature(row)
            feature_names += moa_name
            feature_values += moa_feature
        
        if 'drug_name' in features:
            dr_name, dr_feature = drug_name_feature(row)
            feature_names += dr_name
            feature_values += dr_feature

        if 'monotherapy' in features:
            mono_name, mono_feature = monotherapy_feature(row)
            feature_names += mono_name
            feature_values += mono_feature
        
        if 'chemical_structure' in features:
            ch_name, ch_feature = chemical_structure_feature(row)
            feature_names += ch_name
            feature_values += ch_feature

        if 'molecular' in features:
            include_target_gene = False
            include_network = False
            if 'target_gene' in features:
                include_target_gene = True
                if 'network' in features:
                     include_network = True
            mol_name, mol_feature = molecular_feature(row, include_target_gene, include_network)
            
            feature_names += mol_name
            feature_values += mol_feature

        if 'geneset' in features:
            g_name, g_feature = geneset_feature(row)
            feature_names += g_name
            feature_values += g_feature

        if 'dependency' in features:
            dep_name, dep_feature = cancer_dependency_feature(row)
            feature_names += dep_name
            feature_values += dep_feature

        if 'synthetic_lethality' in features:
            #syn_name, syn_feature = synthetic_lethality_feature(row)
            if not('molecular' in features):
                mol_name, mol_feature = molecular_feature(row, include_target_gene = True, include_network = True)
            syn_name, syn_feature = synthetic_lethality_feature(mol_name, mol_feature)
            feature_names += syn_name
            feature_values += syn_feature    

        return feature_names, feature_values
    
    # initiat feature set(X) and label set(Y)
    X = []
    
    # check if dataset contains label column:
    if return_label and all(re.search('.response_[a-z]+', col) is None for col in data.columns):
            return_label = False
            print('No prediction label target! Return no label.')
    Y = []
    feature_names = []
    # Read by row, construct total feature sets from synergetic experiment information
    data = data.reset_index(drop = True)
    for idx, row in tqdm(data.iterrows(), total = data.shape[0]):
        # feature name list
        feature_names, feature_values = build_features(row, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, features)
        if idx == 0:
            print("Total features:",len(feature_names))

        X.append(feature_values)
        if return_label:
            Y.append(row['.response_'+target])
        
        if if_train:
            switch_cols = {'.metadata_treatment_1':'.metadata_treatment_2',
                    '.metadata_treatment_2':'.metadata_treatment_1',
                    '.metadata_moa_1':'.metadata_moa_2',
                    '.metadata_moa_2':'.metadata_moa_1'}
            new_row = row.rename(switch_cols)
            _, feature_values = build_features(new_row, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, features)
            X.append(feature_values)
            if return_label:
                Y.append(new_row['.response_'+target])
    
    X = np.array(X)
    Y = np.array(Y)

    return feature_names, X, Y



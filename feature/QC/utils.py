#!usr/bin/python3
#author: @rayezh

import pandas as pd
import numpy as np
import collections
from itertools import permutations
import os, json
from scipy import stats

def quantile_normalize_expression_levels(df):
    """
    Parameters:
    df: expression level feature 
    a pandas dataframe

    Yields:
    df_new: quantiled expression level feature
    a pandas dataframe
    """
    cols = df.columns
    exp_cols = [c for c in cols if c.endswith('_exp')]
    
    df_sorted = []
    # sort every row of df_sorted by ascending order
    for idx, r in df.iterrows():
        sorted_row = sorted(r[exp_cols].to_list())
        #print(sorted_row)
        df_sorted.append(sorted_row)
    #print(np.array(df_sorted))
    df_distr = np.mean(np.array(df_sorted), axis = 0)
    # save the distribution: np.savetxt(df_distr, 'exp_distribution.txt')
    #print(df_distr)
    for idx, r in df.iterrows():
        sorted_row = sorted(r[exp_cols].to_list())
        v2i = {sorted_row[i]:i for i in range(len(sorted_row))} #value to index
        for c in exp_cols:
            #print(df.loc[idx,c], df_distr[v2i[df.loc[idx,c]]])
            df.loc[idx,c] = df_distr[v2i[df.loc[idx,c]]]
    return df

def quantile_normalize_RPPA_levels(df):
    """
    Parameters:
    df: expression level feature 
    a pandas dataframe

    Yields:
    df_new: quantiled expression level feature
    a pandas dataframe
    """
    cols = df.columns
    rppa_cols = [c for c in cols if c.endswith('_RPPA')]
    
    df_sorted = []
    # sort every row of df_sorted by ascending order
    for idx, r in df.iterrows():
        sorted_row = sorted(r[rppa_cols].to_list())
        #print(sorted_row)
        df_sorted.append(sorted_row)
    #print(np.array(df_sorted))
    df_distr = np.mean(np.array(df_sorted), axis = 0)
    # save the distribution: np.savetxt(df_distr, 'exp_distribution.txt')
    #print(df_distr)
    for idx, r in df.iterrows():
        sorted_row = sorted(r[rppa_cols].to_list())
        v2i = {sorted_row[i]:i for i in range(len(sorted_row))} #value to index
        for c in rppa_cols:
            #print(df.loc[idx,c], df_distr[v2i[df.loc[idx,c]]])
            df.loc[idx,c] = df_distr[v2i[df.loc[idx,c]]]
    return df


def check_reproducibility(data, datatype, pred_target):
    """ Measure reproducibility of prediction target on datasets
    
    Parameters:
    -----------
    data: Pandas dataframe
        response data (synergy/monotherapy)
    datatype: string
        response data type (synergy/monotherapy)
    
    pred_target: string
        'aoc' or 'bliss'
    Yields:
    -------

    """
    data = data[data['.response_'+pred_target].notna()]

    if datatype == 'synergy':
        data['.metadata_treatment_comb'] = ['_'.join(sorted(list(x))) for _,x in data[['.metadata_treatment_1','.metadata_treatment_2']].iterrows()]
        data['.metadata_treatment_remarks'] = data['.metadata_treatment_remarks'].fillna('')
        experiments = ['_'.join(x) for _,x in data[['.identifier_sample_name',  '.metadata_treatment_comb', '.metadata_treatment_remarks']].iterrows()]
    elif datatype == 'monotherapy':
        data['.metadata_treatment_remarks'] = data['.metadata_treatment_remarks'].fillna('')
        experiments = ['_'.join(x) for _,x in data[['.identifier_sample_name', '.metadata_treatment','.metadata_treatment_remarks']].iterrows()]

    
    data['experiment'] = experiments
    y = collections.Counter(experiments)
    
    # find experiments with replicates
    rep_experiments = [{i:j} for i,j in y.items() if j > 1]
    rep = []
    freq = []
    for i in rep_experiments:
        rep.extend(list(i.keys()))
        freq.extend(list(i.values()))
    
    # check reproducibility
    freq_summary = pd.DataFrame({"experiment": rep, "freq": freq})  #replicates
    print("number of replicates:", freq_summary.shape[0])
    
    score_pair =[]
    for rep in freq_summary['experiment']:
        rows = data[data['experiment'] == rep]
        r,c = rows.shape
        score = list(rows['.response_'+pred_target])
        score_pair.extend([[x,y] for x, y in permutations(score,2)])
        
    score_pair = np.array(score_pair)
    
    try:
        cor_pearson, p_pearson = stats.pearsonr(score_pair[:,0], score_pair[:,1])
        print("Reproducibility (Pearson's cor) of "+pred_target+" across replicated experiments on "+datatype+" dataset is", cor_pearson, ", p =", p_pearson)
        cor_spearman, p_spearman = stats.spearmanr(score_pair[:,0], score_pair[:,1])
        print("Reproducibility (Spearman's cor) of "+pred_target+" across replicated experiments on "+datatype+" dataset is", cor_spearman, ", p =", p_spearman)
    except:
        print("no replicates!")
        cor_pearson, p_pearson = 'nan', 'nan'
        cor_spearman, p_spearman = 'nan', 'nan'
    
    return freq_summary, cor_pearson, p_pearson, cor_spearman, p_spearman
    
def split_monotherapy_synergy(response_df):
    """ Split monotherapy and dual-treatment synergy experimental datasets.
    
    Monotherapy: will be used to construct monotherapy feature sets, used to predict synergy.
    Synergy: prediction target dataset.
    
    Parameters:
    -----------
    response_df: Pandas dataframe

    Yields:
    -------
    cleaned_response_df: Pandas dataframe
    monotherapy_df: Pandas dataframe
    synergy_df: Pandas dataframe
    """
    # Clean out datasets with missing treatment information
    cleaned_response_df = response_df.dropna(subset = ['.metadata_treatment_1', '.metadata_treatment_2'])

    # Select all monotherapy datasets
    """ Monotherapy columns:
        .identifier_batch       
        .identifier_sample_name
        .metadata_moa
        .metadata_treatment_remarks
        .metadata_treatment
        .response_aoc
        .metadata_cancer_type
        .metadata_cancer_subtype
    """
    col_rename = {'.metadata_treatment_1': '.metadata_treatment',
            '.metadata_treatment_2': '.metadata_treatment',
            '.metadata_moa_1' : '.metadata_moa',
            '.metadata_moa_2' : '.metadata_moa'}
    df1 = cleaned_response_df.loc[cleaned_response_df['.metadata_treatment_1'] == "Monotherapy"].drop(columns = ['.metadata_treatment_1', 
        '.metadata_moa_1', 
        '.response_bliss']).rename(columns = col_rename)
    df2 = cleaned_response_df.loc[cleaned_response_df['.metadata_treatment_2'] == "Monotherapy"].drop(columns = ['.metadata_treatment_2', 
        '.metadata_moa_2', 
        '.response_bliss']).rename(columns = col_rename)
    monotherapy_df = pd.concat([df1, df2])
    
    # Select all synergy datasets (for 5-fold cross validation)
    """ Synergy columns:
        .identifier_batch
        .identifier_sample_name
        .metadata_moa_1
        .metadata_moa_2
        .metadata_treatment_remarks
        .metadata_treatment_1
        .metadata_treatment_2
        .response_aoc
        .response_bliss
        .metadata_cancer_type
        .metadata_cancer_subtype
    """
    synergy_df = cleaned_response_df[(cleaned_response_df['.metadata_treatment_1'] != "Monotherapy") & (cleaned_response_df['.metadata_treatment_2'] != "Monotherapy")]

    return cleaned_response_df, monotherapy_df, synergy_df

def construct_monotherapy_response_features(df):
    """ Construct Monotherapy response features
    Monotherapy features: drug responses (.response_aoc) from batch Oncolead_000-012 
    Parameters:
    -----------
    df: Pandas dataframe
        Monotherapy response experimental datasets. 
        Columns:
            .identifier_batch
            .identifier_sample_name
            .metadata_moa
            .metadata_treatment_remarks
            .metadata_treatment
            .response_aoc
            .metadata_cancer_type
            .metadata_cancer_subtype
    
    Yields:
    -------
    df_new: Pandas dataframe
        Columns:
            cell_line.treatment
            response_Oncolead_000
            ...
            response_Oncolead_012

    """
    df['cell_line.treatment'] = df['.identifier_sample_name']+'.'+df['.metadata_treatment']
    all_cell_line_treatment = sorted(set(df['cell_line.treatment']))
    all_batch =  sorted(set(df['.identifier_batch']))
    df_new = {}
    for idx, c_t in enumerate(all_cell_line_treatment):
        r_new = {**{'cell_line.treatment':c_t, **{b:[] for b in all_batch}}}
        responses = df[df['cell_line.treatment'] == c_t].groupby(by = '.identifier_batch')['.response_aoc'].mean()
        for b in all_batch:
            try:
                r_new[b] = responses.loc[b]
            except:
                r_new[b] = np.nan
        df_new.update({idx:r_new})
    df_new = pd.DataFrame.from_dict(df_new, orient = 'index')
    return df_new

def construct_monotherapy_drug_response_network(df):
    """ Construct Monotherapy drug response network
    Monotherapy features: drug responses (.response_aoc) from batch Oncolead_000-012 
    Parameters:
    -----------
    df: Pandas dataframe
        Monotherapy response experimental datasets. 
        Columns:
            .identifier_batch
            .identifier_sample_name
            .metadata_moa
            .metadata_treatment_remarks
            .metadata_treatment
            .response_aoc
            .metadata_cancer_type
            .metadata_cancer_subtype
    
    Yields:
    -------
    drug_network_features

    """
    import networkx as nx 

    all_drugs  = df['.metadata_treatment'].unique()
    # compute median of AoC for each drug-cell line combination
    df_mono_median = df[['.identifier_sample_name', '.metadata_treatment','.response_aoc']].groupby(['.identifier_sample_name', '.metadata_treatment']).agg(np.median).reset_index()
    df_mono_median['.response_aoc'] = df_mono_median.groupby('.identifier_sample_name')['.response_aoc'].transform(lambda x: (x-x.min())/(x.max()-x.min()))
    
    distance_matrix_mono = np.zeros((len(all_drugs), len(all_drugs)))
    removed_edges_mono = []
    distance_matrix_mono_dict = {}
    drug_rwr_features = {}
    drug_rwr_features_dict = {}

    for d1 in all_drugs:
        all_d = []
        for d2 in all_drugs:
            df1 = df_mono_median[df_mono_median['.metadata_treatment'] == d1].rename(columns = {'.response_aoc':'aoc_1', '.metadata_treatment':'treatment_1'})
            df2 = df_mono_median[df_mono_median['.metadata_treatment'] == d2].rename(columns = {'.response_aoc':'aoc_2', '.metadata_treatment':'treatment_2'})
            tmp = df1.merge(df2, left_on = '.identifier_sample_name', right_on = '.identifier_sample_name', how = 'inner').dropna(axis = 0)
            if d1==d2:
                d = 0
                removed_edges_mono.append((d1, d2))
            elif tmp.shape[0] >1:
                d,_ = stats.pearsonr(tmp['aoc_1'], tmp['aoc_2'])
                d = 1-d
            else:
                d = np.NaN
                removed_edges_mono.append((d1, d2))

            distance_matrix_mono_dict[(d1, d2)] = d
            distance_matrix_mono[all_drugs == d1, all_drugs == d2] = d
            all_d.append(d)
        drug_rwr_features.update({d1:all_d})

    G = nx.DiGraph(np.array(distance_matrix_mono))
    # relabel the nodes
    mapping = {i:all_drugs[i] for i in range(len(all_drugs))}
    G = nx.relabel_nodes(G, mapping)
    for (u,v) in removed_edges_mono:
        try:
            G.remove_edge(u, v)
        except:
            pass

    # use random walk with restart to calculate the embeddings 
    alpha = 0.05
    aff = []
    for drug in all_drugs:
        drug_rwr_features.update({drug:list(nx.pagerank(G, alpha, personalization={drug:1}, weight = 'weight').values())})
        drug_rwr_features_dict.update({drug:dict(nx.pagerank(G, alpha, personalization={drug:1}, weight = 'weight'))})

    return drug_rwr_features, drug_rwr_features_dict, all_drugs


def map_cell_line():
    """ Mapping cell lines' DepMap ID to ONCOLEAD name

    Parameters:
    -----------
    Yields:
    -------
    cell2id: a dict
        {DepMap_ID:ONCOLEAD_CELL_<name>}
    """
    # Load DepMap cell line mapping information
    cell_map = pd.read_csv('../cancer_dependency_features/sample_info.csv', header = 0)

    # Load all cell lines in DDR drug synergy experiments
    data = pd.read_csv('../data/synergy_responses.tsv', sep = '\t', header = 0)
    all_cell = list(set(data['.identifier_sample_name']))

    # build cell line mapping
    id2cell = {}  #dictionary: ONCOLEAD_cell:DepMap_ID

    for c in all_cell:
        name = c.split('_')[-1]
        if name in list(cell_map['stripped_cell_line_name']):
            cid = cell_map.loc[cell_map['stripped_cell_line_name'] == name, 'DepMap_ID'].iloc[0]
            id2cell.update({cid:c})
        else:
            print("Missing cell line in DepMap:", name)

    return id2cell

def build_moa_cancer_dependency(dep):
    """ Get cancer dependency of drugs based on target genes

    Parameters:
    -----------
    dep: Pandas dataframe;
        columns: DepMap_ID of cell lines
        index: gene name (Hugo ID)
    Yields:
    -------
    all_moa_dep: dictionary
        key1: ONCOLEAD cell name
        key2: moa (mode of action)
        key3: max/mean/min
    """

    drug2gene = json.load(open('../target_gene/drug2gene.json'))  # check from here
    id2cell = map_cell_line()

    all_moa_dep = {}
    for moa, genes in drug2gene.items():
        moa_dep = {}
        for dep_id in dep.columns:
            if dep_id in id2cell.keys():
                cell = id2cell[dep_id]
                gene_dep = []
                for g in genes:
                    try:
                        d = dep.loc[g,dep_id]
                        gene_dep.append(d)
                    except:
                        print(g)
                
                if len(gene_dep)>0:
                    all_dep = {'max':max(gene_dep), 'mean': np.mean(gene_dep), 'min':min(gene_dep)} #select max dependency of all target genes
                else:
                    all_dep = {'max':np.nan, 'mean':np.nan, 'min':np.nan}
                
                moa_dep.update({cell:all_dep})
            else:
                pass
        all_moa_dep.update({moa:moa_dep})

    return all_moa_dep



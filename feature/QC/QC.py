#!usr/bin/python3
import pandas as pd
import numpy as np
import collections
from collections import defaultdict
import os, json, pickle
from tqdm import tqdm
from scipy import stats
from glob import glob
from itertools import combinations
from venn import venn
from matplotlib import pyplot 
from utils import *

# chemical structure required modules
from rdkit import Chem,DataStructs
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdmolops
from pubchempy import *
import openbabel
import pybel

def geneset_qc():
    """ QC on geneset information; Construct geneset feature matrix 
    
    The geneset matrix will be laid out like:
         |geneset1|geneset2|geneset3|...
    -------------------------------------
    gene1|   0    |    1   |   1    |
    -------------------------------------
    gene2|   0    |    0   |   1    |
    -------------------------------------
    gene3|   1    |    1   |   1    |
    -------------------------------------
    ...  |

    Params
    -------

    Yields
    ------
    """
    df = pd.read_csv('../data/genesets.tsv', sep = '\t', header =0, index_col = 0)
    print("QC on geneset information ...")
    
    outpath = '../geneset_features/'
    os.makedirs(outpath, exist_ok = True)
    
    #print geneset total
    geneset_all = list(set(df['geneset_name'].to_list()))
    print("Total genesets:",len(geneset_all)) #9450
    f = open(outpath+'all_geneset.txt','w')
    f.write('\n'.join(geneset_all))
    f.close()

    #print gene total
    gene_all = sorted(list(set(df['gene_symbol'].to_list())))
    print("Total genes:",len(gene_all)) #18681
    df_mol = pd.read_csv('../molecular_features/exp_molecular_features.tsv', sep = '\t', header =0, index_col = 0)
    gene = [x.split('_')[0] for x in df_mol.columns[3:]]
    overlap_gene = list(set(gene_all) &set(gene))

    # binray feature: 0 or 1
    features_df = {geneset:[] for geneset in geneset_all}
    for geneset in tqdm(geneset_all):
        gene_geneset = df.loc[df['geneset_name'] == geneset]['gene_symbol'].to_list()
        features_df[geneset] = [1 if x in gene_geneset else 0 for x in overlap_gene]
    features_df = pd.DataFrame.from_dict(features_df, orient = 'index', columns = overlap_gene)
    print(features_df)
    features_df.to_csv(outpath+"geneset_features.csv")

def molecular_qc():
    """
    Data types:

    1. Single gene molecular features: <GENESYMBOL>_snv/cnv/exp/lof;
    2. Geneset features: <GENESET>__coh_pat/lof_pat (geneset sep by '_');
    3. DDR features: *_ddr (78);
    
    Params
    ------
    Yields
    ------
    """

    # Read molecular raw data
    path = '../data/molecular_features.tsv'
    data = pd.read_csv(path,sep = '\t',index_col=0 ,header = 0)
    
    # Create output path
    out_path = '../molecular_features/'
    print("Saving processed molecular features files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)

    # Quantile normalize gene expression level
    print("Quantile normalize gene expression levels of molecular features ...")
    data = quantile_normalize_expression_levels(data)
    data = data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    data.to_csv(out_path+'quantiled_molecular_features.tsv', index = False, sep = '\t')

    
    # Split datasets
    print("Split molecular data sets ...")
    columns = list(data.columns)
    snv = []
    cnv = []
    exp = []
    lof = []

    coh_pat = []
    lof_pat = []

    ddr = []
    uk = []
    while len(columns)>0:
        col = columns.pop()
        if col.endswith('_snv'):
            snv.append(col)
        elif col.endswith('_cnv'):
            cnv.append(col)
        elif col.endswith('_exp'):
            exp.append(col)
        elif col.endswith('_lof'):
            lof.append(col)
        elif col.endswith('__coh_pat'):
            coh_pat.append(col)
        elif col.endswith('__lof_pat'):
            lof_pat.append(col)
        elif col.endswith('_ddr'):
            ddr.append(col)
        else:
            uk.append(col)

    print("snv:", len(snv))
    print("cnv:", len(cnv))
    print("exp:", len(exp))
    print("lof:", len(lof))
    print("coh_pat:", len(coh_pat))
    print("lof_pat:", len(lof_pat))
    print("ddr:", len(ddr))
    print("uk columns:", len(uk), uk)
    celltype = ['.identifier_sample_name']
    
    print("intersected genes between exp, snv, snv and lof: to mol.pdf")
    all_mols = {"snv":set([x.split("_")[0] for x in snv]), "cnv":set([x.split("_")[0] for x in cnv]), "exp":set([x.split("_")[0] for x in exp]), "lof":set([x.split("_")[0] for x in lof])}
    #print(all_mols)
    plt = venn(all_mols)
    pyplot.savefig("mol.pdf")

    #exp: QC of gene expression level
    exp_cols = celltype+exp
    exp_data = data[exp_cols]
    #exp_data = exp_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    exp_data.to_csv(out_path+'exp_molecular_features.tsv', index = False, sep = '\t')

    # make dictionary: gene, exp -> rank of exp in thia gene
    dic={}
    for col in exp:
        gene = col.split('_')[0]
        md = np.median(exp_data[col].to_list())
        mn = np.mean(exp_data[col].to_list())
        sd =  np.std(exp_data[col].to_list())
        d = {'median':md, 'mean': mn,'sd': sd}
        #values = sorted(set(exp_data[col].to_list()))
        #d = {np.round(v, 4):i for i, v in enumerate(values)}
        dic.update({gene:d})
    print(dic)

    # TODO
    pickle.dump(dic, open('exp_rank.pkl', 'wb'))
    # TODO
    ddr_gene_f = open('../target_gene/ddr_gene.txt', 'w') 
    ddr_gene_f.write('\n'.join(dic.keys()))
    ddr_gene_f.close()

    #snv: QC of number of alleles affected by pathigenic mutation or indel
    snv_cols = celltype+snv
    snv_data = data[snv_cols]
    #snv_data = snv_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    snv_data.to_csv(out_path+'snv_molecular_features.tsv', index = False, sep = '\t')
    
    #cnv: QC of gene copy number variation
    cnv_cols = celltype+cnv
    cnv_data = data[cnv_cols]
    #cnv_data = cnv_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    cnv_data.to_csv(out_path+'cnv_molecular_features.tsv', index = False, sep = '\t')

    #lof: QC of predicted gene loss of function
    lof_cols = celltype+lof
    lof_data = data[lof_cols]
    #lof_data = lof_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    lof_data.to_csv(out_path+'lof_molecular_features.tsv', index = False, sep = '\t')

    # coh_pat
    coh_pat_cols = celltype+coh_pat
    coh_pat_data = data[coh_pat_cols]
    #coh_pat_data = coh_pat_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    coh_pat_data.to_csv(out_path+'coh_pat_molecular_features.tsv', index = False, sep = '\t')
    
    # lof_pat
    lof_pat_cols = celltype+lof_pat
    lof_pat_data = data[lof_pat_cols]
    #lof_pat_data = lof_pat_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    lof_pat_data.to_csv(out_path+'lof_pat_molecular_features.tsv', index = False, sep = '\t')
    
    # ddr
    ddr_cols = celltype + ddr
    ddr_data = data[ddr_cols]
    #ddr_data = ddr_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    ddr_data.to_csv(out_path+'ddr_molecular_features.tsv', index = False, sep = '\t')

    # RPPA: TODO
    #data = pd.read_csv('../molecular_features/RPPA.tsv', sep = '\t')
    data = quantile_normalize_expression_levels(data)


def response_qc(path, fname):
    """ QC of drug response data

    * publication 1 version: 
        1. use synergy that must with monotherapy response data
        2. remove 'MSC*' in drug names

    Params
    ------
    Yields
    ------
    """
    def clean_missing_monotherapy(mono_features, synergy):
        """ Remove synergy responses with missing monotherapy responses data
        
        Params
        ------
        mono_features: Pandas dataframe
        synergy: Pandas dataframe
        
        Yields
        ------
        synergy: synergy after removed samples without monotherapy response features
        """

        idx = []
        for i, r in synergy.iterrows():
            if ((r['.identifier_sample_name']+'.'+r['.metadata_treatment_1'] in list(mono_features['cell_line.treatment'])) and (r['.identifier_sample_name']+'.'+r['.metadata_treatment_2'] in list(mono_features['cell_line.treatment']))):
                idx.append(i)
        synergy = synergy.loc[idx]
        return synergy
    
    
    
    # preprocess: remove used drugs
    path = '../data/'
    df = pd.read_csv(path+fname, sep = '\t', header = 0,index_col = 0)
    print(df)

    # drop nan
    df = df.dropna(subset=['.metadata_treatment_1', '.metadata_treatment_2'])
    
    # preprocess: remove 'MSC**' in drug names
    df[['.metadata_treatment_1', '.metadata_treatment_2']] = df[['.metadata_treatment_1', '.metadata_treatment_2']].replace('MSC[A-Z0-9]+', '', regex=True).replace('_+', '_', regex=True).replace('^_', '', regex=True).replace('_$', '', regex=True)
    
    # replace moa subclasses with parent class
    responses =  df.replace('Cytostatic_Antimetabolite_Pyrimidineanalogue', 'Cytostatic_Antimetabolite').replace('Cytostatic_Alkylator', 'Cytostatic_Intercalator').replace('MEKi', 'CDKi').replace('CDK4i_6i', 'CDKi')
    
    # split monotherapy and synergy data
    responses, monotherapy, synergy = split_monotherapy_synergy(responses)
    
    # save cleaned/splited datasets
    responses.to_csv(path+'cleaned_'+fname, sep = '\t', index = False,  na_rep = 'NA')
    monotherapy.to_csv(path+'monotherapy_responses.tsv', sep ='\t', index = False,  na_rep = 'NA')
    synergy.to_csv(path+'synergy_responses.tsv', sep ='\t', index = False,  na_rep = 'NA')

    # construct monotherapy response features
    out_path = '../monotherapy_features/'
    print("Saving processed monotherapy features files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)
    mono_features = construct_monotherapy_response_features(monotherapy)
    mono_features.to_csv(out_path+'monotherapy_features.tsv', sep = '\t', index = False, na_rep = 'NA')

    # clean synergy(remove samples missing monotherapy responses)
    print("Before cleaning:", synergy.shape[0])
    synergy = clean_missing_monotherapy(mono_features, synergy)
    print("After cleaning:", synergy.shape[0])
    synergy.to_csv(path+'synergy_responses_with_monotherapy.tsv', sep ='\t', index = False,  na_rep = 'NA')

    # check #treatment in Monotherapy
    treat = list(set(monotherapy['.metadata_treatment']))
    print("Total treatment in monotherapy:", len(treat))
    #print(treat)

    # check #treatment in dual-therapy
    treat = list(set(list(synergy['.metadata_treatment_1'])+list(synergy['.metadata_treatment_2'])))
    print("Total treatment in dual-therapy:", len(treat))
    #print(treat)
    
    # check #cell line
    cl = list(set(monotherapy['.identifier_sample_name']))
    print("Total cell lines in monotherapy:", len(cl))
    cl = list(set(synergy['.identifier_sample_name']))
    print("Total cell lines in synergy:", len(cl))
    
    # check #tissue   
    tissue = list(set(monotherapy['.metadata_cancer_type']))
    print("Total cell lines in monotherapy:", len(tissue))
    tissue = list(set(synergy['.metadata_cancer_type']))
    print("Total cell lines in synergy:", len(tissue))

    # check overall preprodicibility
    os.makedirs('./reproducibility', exist_ok = True)
    freq_sum = check_reproducibility(responses, 'synergy', 'aoc')
    freq_sum.to_csv("./reproducibility/response_replictats_aoc.tsv", sep ='\t', index = False )
    freq_sum = check_reproducibility(responses, 'synergy', 'bliss')
    freq_sum.to_csv("./reproducibility/response_replictats_bliss.tsv", sep ='\t', index = False )
       
    # check synergy reprodicibility
    freq_sum = check_reproducibility(synergy, 'synergy', 'aoc')
    freq_sum.to_csv("./reproducibility/synergy_replictats_aoc.tsv", sep ='\t', index = False )
    freq_sum = check_reproducibility(synergy, 'synergy', 'bliss')
    freq_sum.to_csv("./reproducibility/synergy_replictats_bliss.tsv", sep ='\t', index = False )

    # check in-batch:
    for batch in sorted(set(synergy['.identifier_batch'])):
        print(batch)
        synergy_new = synergy[synergy['.identifier_batch'] == batch]
        _ = check_reproducibility(synergy_new, 'synergy', 'aoc')
        _ = check_reproducibility(synergy_new, 'synergy', 'bliss')

    # check monotherapy reprodicibility
    freq_sum = check_reproducibility(monotherapy, 'monotherapy', 'aoc')
    freq_sum.to_csv("./reproducibility/monotherapy_replictats_aoc.tsv", sep ='\t', index = False )

    # check in-batch
    for batch in sorted(set(monotherapy['.identifier_batch'])):
        print(batch)
        monotherapy_new = monotherapy[monotherapy['.identifier_batch'] == batch]
        _ = check_reproducibility(monotherapy_new, 'monotherapy', 'aoc')

    # check drug_pair tested in cell lines
    synergy['drug_pair'] = ['.'.join(sorted(list(r))) for _,r in synergy[['.metadata_treatment_1','.metadata_treatment_2']].iterrows()]
    cellline_summary = synergy.groupby('drug_pair')['.identifier_sample_name'].nunique().sort_values()
    mean = cellline_summary.mean()
    print("Averge # cellines tested for drug pairs:", mean)
    print(cellline_summary)

    # check combinations
    synergy['combination'] = synergy['drug_pair']+'.'+synergy['.identifier_sample_name']
    combination_summary = synergy.groupby('combination')['combination'].count().sort_values()
    print(combination_summary)

    # get all drugs tested
    drugs_summary = monotherapy[['.metadata_treatment', '.metadata_moa']].drop_duplicates().sort_values(by = '.metadata_treatment').rename(columns = {'.metadata_treatment':"drug_name",  '.metadata_moa':"mode-of-action"})
    os.makedirs('../target_gene', exist_ok = True)
    drugs_summary.to_csv("../target_gene/all_drugs_summary.csv", index = False)

def encode_categorical_feature():
    """ Categorical feature encoding object
	
	Params
	------
	data: '../data/synergy_responses_with_monotherapy.tsv' pandas frame. the total response dataset
	
	Yields
	------
	batch: dict
	cell_line: dict
	moa: dict
	remark: dict
	cancer_type: dict
	cancer_subtype: dict
	"""
    data = pd.read_csv('../data/synergy_responses_with_monotherapy.tsv', sep = '\t')
	
    outpath = '../categorical_embedding_features/'
    os.makedirs(outpath, exist_ok = True)
    
    batch = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(data['.identifier_batch'])))})
    cell_line = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(data['.identifier_sample_name'])))})
    moa = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set([str(x) for x in data['.metadata_moa_1'].to_list()+data['.metadata_moa_2'].to_list()])))})
    treatment = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set([str(x) for x in data['.metadata_treatment_1'].to_list()+data['.metadata_treatment_2'].to_list()])))})
    remark = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set([str(x) for x in data['.metadata_treatment_remarks']])))})
    cancer_type = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(data['.metadata_cancer_type'])))})
    cancer_subtype = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(data['.metadata_cancer_subtype'])))})
    
    json.dump(batch, open(outpath+'batch.json', 'w'))
    json.dump(cell_line, open(outpath+'cell_line.json', 'w'))
    json.dump(moa, open(outpath+'moa.json', 'w'))
    json.dump(treatment, open(outpath+'treatment.json', 'w'))
    json.dump(remark, open(outpath+'remark.json', 'w'))
    json.dump(cancer_type, open(outpath+'cancer_type.json', 'w'))
    json.dump(cancer_subtype, open(outpath+'cancer_subtype.json', 'w'))

def chemical_structure_qc():
    """ QC of chemical structure information and feature preprocess

    Transform SMILE formatted chemical structure into 6 types of Fingerprints
    """

    def fp_to_feature(fp, max_bit):
        """
        """
        feature = [0 for i in range(max_bit)]
        for n in fp.bits:
            feature[n] = 1
        return feature

    outpath = '../chemical_structure_features/'
    os.makedirs(outpath, exist_ok = True)

    df_map = pd.read_csv('../data/chemical_structures.csv', header = 0)
    df_map = df_map.drop_duplicates(subset = ['treatment_name_orig'])

    # make new dirs
    os.makedirs(outpath+"MACCS_features/", exist_ok = True)
    os.makedirs(outpath+"Morgan_features/", exist_ok = True)
    os.makedirs(outpath+"RDK_features/", exist_ok = True)
    os.makedirs(outpath+"FP2_features/", exist_ok = True)
    os.makedirs(outpath+"FP3_features/", exist_ok = True)
    os.makedirs(outpath+"FP4_features/", exist_ok = True)

    max_bit_fp2 = 1024
    max_bit_fp3 = 56
    max_bit_fp4 = 308

    all_features = []
    for _,r in df_map.iterrows():
        smile = r['treatment_smiles']
        ms = Chem.MolFromSmiles(smile)
        mol = pybel.readstring("smi", smile)
        
        # MACCS features (167*1)
        fp = MACCSkeys.GenMACCSKeys(ms)
        tmp = fp.ToBitString()
        feature_1 = list(map(int, tmp))
        np.savetxt(outpath+'MACCS_features/'+r['treatment_name_orig'],np.array(feature_1))
        feature_1 = {'MACCS_'+str(idx+1):x for idx,x in enumerate(feature_1)}

        # Morgan Fingerprints (1024*1)
        fp = AllChem.GetMorganFingerprintAsBitVect(ms,2,nBits=1024)
        tmp = fp.ToBitString()
        feature_2 = list(map(int, tmp))
        np.savetxt(outpath+'Morgan_features/'+r['treatment_name_orig'],np.array(feature_2))
        feature_2 = {'Morgan_'+str(idx+1):x for idx,x in enumerate(feature_2)}
        
        # FP2 (1024*1)
        fp = mol.calcfp('FP2')
        feature_3 = fp_to_feature(fp, max_bit_fp2)
        np.savetxt(outpath+'FP2_features/'+r['treatment_name_orig'],np.array(feature_3))
        feature_3 = {'FP2_'+str(idx+1):x for idx,x in enumerate(feature_3)}

        # FP3 (56*1)
        fp = mol.calcfp('FP3')
        feature_4 = fp_to_feature(fp, max_bit_fp3)
        np.savetxt(outpath+'FP3_features/'+r['treatment_name_orig'],np.array(feature_4))
        feature_4 = {'FP3_'+str(idx+1):x for idx,x in enumerate(feature_4)}

        # FP4 (308*1)
        fp = mol.calcfp('FP4')
        feature_5 = fp_to_feature(fp, max_bit_fp4)
        np.savetxt(outpath+'FP4_features/'+r['treatment_name_orig'],np.array(feature_5))
        feature_5 = {'FP4_'+str(idx+1):x for idx,x in enumerate(feature_5)}

        # RDK Fingerprints (2048*1)
        fp = rdmolops.RDKFingerprint(ms)
        tmp = fp.ToBitString()
        feature_6 = list(map(int, tmp))
        np.savetxt(outpath+'RDK_features/'+r['treatment_name_orig'],np.array(feature_6))
        feature_6 = {'RDK_'+str(idx+1):x for idx,x in enumerate(feature_6)}

        drug_name = {'treatment':r['treatment_name_orig']}
        the_features = {**drug_name,**feature_1,**feature_2,**feature_3,**feature_4,**feature_5,**feature_6}
        all_features.append(the_features)

    all_feature=pd.DataFrame.from_dict(all_features)
    all_feature.to_csv(outpath+'all_chemical_structure.csv', index = False)

def cancer_dependency_qc():
    """ QC and feature construction of cancer dependency information
    """
    out_path = '../cancer_dependency_features/'
    print("Saving processed cancer dependency feature files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)
    
    # CERES dataset
    dep = pd.read_csv(out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CERES_FC.txt', sep = '\t', index_col = 0)
    moa_dep = build_moa_cancer_dependency(dep)
    p_out = out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CERES_FC_dep.json'
    json.dump(moa_dep, open(p_out, 'w'))

    # CRISPRcleanR datase
    dep = pd.read_csv(out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CRISPRcleanR_FC.txt', sep = '\t', index_col = 0)
    moa_dep = build_moa_cancer_dependency(dep)
    p_out = out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CRISPRcleanR_FC_dep.json'
    json.dump(moa_dep, open(p_out, 'w'))

def main():
    response_qc(path = '../data/', fname = 'responses_training.tsv')
    encode_categorical_feature()
    molecular_qc()
    geneset_qc()
    chemical_structure_qc()
    #cancer_dependency_qc()

if __name__ == "__main__":
    main()

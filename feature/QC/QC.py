#!usr/bin/python3
"""
This program has functions as below:
1. QC on all datasets used in this study
2. Construct feature matrix for each data type
"""

import pandas as pd
import numpy as np
import collections
import os, json, pickle
from tqdm import tqdm
from scipy import stats
from glob import glob
from itertools import combinations
#from venn import venn
from matplotlib import pyplot 
from utils import *

# change drug names
updated_drug = {
'ABT737':'ABT737','ABT888_Veliparib':'veliparib','AG014699_Rucaparib':'rucaparib','AZD4547':'AZD4547','Arabinofuranoside_Cytarabine':'cytarabine','Azacytidine':'azacytidine',
'BI2536':'BI2536','BML284':'BML284','BMN673_Talazoparib':'talazoparib','BMXIN1':'BMXIN1','Bleomycin':'bleomycin','Bortezomib':'bortezomib','CPI0610':'CPI0610',
'Carboplatin':'carboplatin','Carfilzomib':'carfilzomib','Ceritinib':'ceritinib','Cisplatin':'cisplatin','Crizotinib':'crizotinib','Dacarbazine':'dacarbazine',
'Doxorubicin':'doxorubicin','Duvelisib':'duvelisib','Erlotinib':'erlotinib','Etoposide':'etoposide','Everolimus':'everolimus','Fluorouracil':'fluorouracil','Formestane':'formestane',
'GAMMA':'gamma','GDC0425':'GDC0425','GDC0941_Pictilisib':'pictilisib','GS9973_Entospletinib':'entospletinib','GSK2334470':'GSK2334470','GSK429286A':'GSK429286A','GSK525762A':'GSK525762A',
'GSK923295':'GSK923295','Geldanamycin':'geldanamycin','Gemcitabine':'gemcitabine','Gilteritinib':'gilteritinib','IBET151':'IBET151','Ibrutinib':'ibrutinib','Imatinib':'imatinib',
'K858':'K858','LY2228820_Ralimetinib':'ralimetinib','LY2584702':'LY2584702','LY2603618':'LY2603618','LY2606368_Prexasertib':'prexasertib','LY2857785':'LY2857785','Lapatinib':'lapatinib',
'Lomustine':'lomustine','M3541':'M3541','M3814':'peposertib','M4076':'M4076','MK4827_Niraparib':'niraparib','ML167':'ML167', 'Methotrexate': 'methotrexate','Melphalan':'melphalan','Mitomycin':'mitomycin','NSC348884':'NSC348884',
'Nutlin':'nutlin','OSI027':'OSI027','Olaparib':'olaparib','Oxaliplatin':'oxaliplatin','PD0332991_Palbociclib':'palbociclib','PF04217903':'PF04217903','PLX647':'PLX647','Paclitaxel':'paclitaxel',
'Pomalidomide':'pomalidomide','Purvalanol':'purvalanol','RO3306':'RO3306','Roscovitine_Seliciclib':'seliciclib','Ruxolitinib':'ruxolitinib','SCH727965_MK7965_Dinaciclib':'dinaciclib','SGX523':'SGX523',
'SN38':'SN38','Sorafenib':'sorafenib','Sunitinib':'sunitinib','THZ1':'THZ1','Temozolomide':'temozolomide','Thiotepa':'thiotepa','Topotecan':'topotecan','Ulixertinib':'ulixertinib','VX400_M1774':'M1774',
'VX680_Tozasertib':'tozasertib','VX803_M4344':'gartisertib','VX970_Berzosertib_M6620':'berzosertib','Vatalanib':'vatalanib','Vemurafenib':'vemurafenib','Vincristine':'vincristine','XL184_Cabozantinib':'cabozantinib'
}



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
    df = pd.read_csv('../data/merck_confidential_genesets_20200414.tsv', sep = '\t', header =0, index_col = 0)
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
    df_mol = pd.read_csv('../data/exp_merck_confidential_features_20200221.tsv', sep = '\t', header =0, index_col = 0)
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
    path = '../data/merck_confidential_features_20200221.tsv'
    data = pd.read_csv(path,sep = '\t',index_col=0 ,header = 0)
    data['.identifier_sample_name'] = data['.identifier_sample_name'].apply(lambda x: x.replace('ONCOLEAD_CELL_','')) # remove 'ONCOLEAD_CELL_' in cell line names

    
    # Create output path
    out_path = '../molecular_features/'
    print("Saving processed molecular features files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)

    # Quantile normalize gene expression level
    print("Quantile normalize gene expression levels of molecular features ...")
    data = quantile_normalize_expression_levels(data)
    data = data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    data.to_csv(out_path+'quantiled_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

    
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
    all_mols = {"snv":set([x.split("_")[0] for x in snv]), "cnv":set([x.split("_")[0] for x in cnv]), "exp":set([x.split("_")[0] for x in exp]), "lof":set([x.split("_")[0] for x in lof])}
    #print(all_mols)
    plt = venn(all_mols)
    pyplot.savefig("mol.pdf")

    #exp: QC of gene expression level
    exp_cols = celltype+exp
    exp_data = data[exp_cols]
    #exp_data = exp_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    exp_data.to_csv(out_path+'exp_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

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
    #print(dic)
    pickle.dump(dic, open('exp_rank.pkl', 'wb'))
    ddr_gene_f = open('ddr_gene.txt', 'w') 
    ddr_gene_f.write('\n'.join(dic.keys()))
    ddr_gene_f.close()

    #snv: QC of number of alleles affected by pathigenic mutation or indel
    snv_cols = celltype+snv
    snv_data = data[snv_cols]
    #snv_data = snv_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    snv_data.to_csv(out_path+'snv_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
    
    #cnv: QC of gene copy number variation
    cnv_cols = celltype+cnv
    cnv_data = data[cnv_cols]
    #cnv_data = cnv_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    cnv_data.to_csv(out_path+'cnv_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

    #lof: QC of predicted gene loss of function
    lof_cols = celltype+lof
    lof_data = data[lof_cols]
    #lof_data = lof_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    lof_data.to_csv(out_path+'lof_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

    # coh_pat
    coh_pat_cols = celltype+coh_pat
    coh_pat_data = data[coh_pat_cols]
    #coh_pat_data = coh_pat_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    coh_pat_data.to_csv(out_path+'coh_pat_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
    
    # lof_pat
    lof_pat_cols = celltype+lof_pat
    lof_pat_data = data[lof_pat_cols]
    #lof_pat_data = lof_pat_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    lof_pat_data.to_csv(out_path+'lof_pat_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
    
    # ddr
    ddr_cols = celltype + ddr
    ddr_data = data[ddr_cols]
    #ddr_data = ddr_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    ddr_data.to_csv(out_path+'ddr_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

    # RPPA
    # TODO: add this to all molecular features
    data = pd.read_csv('../molecular_features/RPPA.tsv', sep = '\t')
    data['.identifier_sample_name'] = data['.identifier_sample_name'].apply(lambda x: x.replace('ONCOLEAD_CELL_','')) # remove 'ONCOLEAD_CELL_' in cell line names
    data = quantile_normalize_RPPA_levels(data)
    data.to_csv(out_path+'quantiled_rppa_merck_confidential_features_20200221.tsv', index = False, sep = '\t')


def response_qc(path, fname):
    """ QC of drug response data

    * publication 1 version: 
        1. use synergy that must with monotherapy response data
        2. remove 'MSC*' in drug names
        3. remove th following compounds:
            M2698
            M7583
            Tepotinib
            M8891
            SRA737
            MSC2531331
            AZD6738
            BAY1895344
    
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
        1. cleaned responses
        2. monotherapy responses
        3. synergy: synergy after removed samples without monotherapy response features
        4. reproducibility check in the dataset
        """

        idx = []
        for i, r in synergy.iterrows():
            if ((r['.identifier_sample_name']+'.'+r['.metadata_treatment_1'] in list(mono_features['cell_line.treatment'])) and (r['.identifier_sample_name']+'.'+r['.metadata_treatment_2'] in list(mono_features['cell_line.treatment']))):
                idx.append(i)
        synergy = synergy.loc[idx]
        return synergy
    
    
    # preprocess: remove used drugs
    df = pd.read_csv(path+fname, sep = '\t', header = 0,index_col = 0)
    print(df)

    # drop nan
    df = df.dropna(subset=['.metadata_treatment_1', '.metadata_treatment_2'])
    # exclude drugs
    df = df.loc[(df['.metadata_treatment_1'] != 'MSC2531331')&(df['.metadata_treatment_1'] != 'MSC2691088_SRA737')&(df['.metadata_treatment_2'] != 'MSC2531331')&(df['.metadata_treatment_2'] != 'MSC2691088_SRA737')].copy().reset_index(drop = True)
    # preprocess: remove 'MSC**' in drug names
    df[['.metadata_treatment_1', '.metadata_treatment_2']] = df[['.metadata_treatment_1', '.metadata_treatment_2']].replace('MSC[A-Z0-9]+', '', regex=True).replace('_+', '_', regex=True).replace('^_', '', regex=True).replace('_$', '', regex=True)
    # preprocess: remove the following drugs: M2689, M7583, Tepotinib, M8891
    df = df.loc[~df['.metadata_treatment_1'].str.contains('M2698|M7583|Tepotinib|M8891|AZD6738|BAY1895344')]
    df = df.loc[~df['.metadata_treatment_2'].str.contains('M2698|M7583|Tepotinib|M8891|AZD6738|BAY1895344')]
    
    # replace moa subclasses with parent class
    responses =  df.replace('Cytostatic_Antimetabolite_Pyrimidineanalogue', 'Cytostatic_Antimetabolite').replace('Cytostatic_Alkylator', 'Cytostatic_Intercalator').replace('MEKi', 'CDKi').replace('CDK4i_6i', 'CDKi')
    responses = responses.replace(updated_drug) # replace drug names with updated names
    responses['.identifier_sample_name'] = responses['.identifier_sample_name'].apply(lambda x: x.replace('ONCOLEAD_CELL_','')) # remove 'ONCOLEAD_CELL_' in cell line names
    
    # save processed responses
    responses.to_csv(path+'processed'+fname, sep = '\t')

    # split monotherapy and synergy data
    responses, monotherapy, synergy = split_monotherapy_synergy(responses)
    
    # save cleaned/splited datasets
    responses.to_csv(path+'cleaned_'+fname, sep = '\t', index = False,  na_rep = 'NA') 
    monotherapy.to_csv(path+'monotherapy_responses.tsv', sep ='\t', index = False,  na_rep = 'NA')
    synergy.to_csv(path+'synergy_responses.tsv', sep ='\t', index = False,  na_rep = 'NA') # all synergy responses; could contain samples without monotherapy response features

    # construct monotherapy response features
    out_path = '../monotherapy_features/'
    print("Saving processed monotherapy features files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)
    mono_features = construct_monotherapy_response_features(monotherapy)
    if 'ho1' in fname:
        mono_features.to_csv(out_path+'monotherapy_features_ho1.tsv', sep = '\t', index = False, na_rep = 'NA')
    elif 'ho2' in fname:
        mono_features.to_csv(out_path+'monotherapy_features_ho2.tsv', sep = '\t', index = False, na_rep = 'NA')
    else:
        mono_features.to_csv(out_path+'monotherapy_features.tsv', sep = '\t', index = False, na_rep = 'NA')

    # Construct drug-drug similarity rwr features
    # TODO: compare this feature with other drug features- chemical similarity, drug target similarity, etc.
    # this feature can be generated data-wise
    out_path = '../drug_similarity_features/'
    print("Saving processed drug similarity features files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)
    drug_rwr_features, drug_rwr_features_dict, all_drugs  = construct_monotherapy_drug_response_network(monotherapy)
    if 'ho1' in fname:
        pickle.dump(drug_rwr_features, open(out_path+'drug_similarity_features_ho1.pkl', 'wb'))
        pickle.dump(drug_rwr_features_dict, open(out_path+'drug_similarity_features_dict_ho1.pkl', 'wb'))
        # save all drugs to a file with index
        all_drugs_tmp = pd.DataFrame.from_dict({'drug': all_drugs, 'index': list(range(len(all_drugs)))})
        all_drugs_tmp.to_csv(out_path+'all_drugs_ho1.tsv', sep = '\t', index = True, na_rep = 'NA')
    elif 'ho2' in fname:
        pickle.dump(drug_rwr_features, open(out_path+'drug_similarity_features_ho2.pkl', 'wb'))
        pickle.dump(drug_rwr_features_dict, open(out_path+'drug_similarity_features_dict_ho2.pkl', 'wb'))
        all_drugs_tmp = pd.DataFrame.from_dict({'drug': all_drugs, 'index': list(range(len(all_drugs)))})
        all_drugs_tmp.to_csv(out_path+'all_drugs_ho1.tsv', sep = '\t', index = True, na_rep = 'NA')
    else:
        pickle.dump(drug_rwr_features, open(out_path+'drug_similarity_features.pkl', 'wb'))
        pickle.dump(drug_rwr_features_dict, open(out_path+'drug_similarity_features_dict.pkl', 'wb'))
        all_drugs_tmp = pd.DataFrame.from_dict({'drug': all_drugs, 'index': list(range(len(all_drugs)))})
        all_drugs_tmp.to_csv(out_path+'all_drugs.tsv', sep = '\t', index = True, na_rep = 'NA')    

    # clean synergy(remove samples missing monotherapy responses)
    print("Before cleaning:", synergy.shape[0])
    synergy = clean_missing_monotherapy(mono_features, synergy)
    print("After cleaning:", synergy.shape[0])
    synergy.to_csv(path+'synergy_responses_with_monotherapy.tsv', sep ='\t', index = False,  na_rep = 'NA') # synergy responses with monoterapy response features; used for training models

    # QC reports
    rep_path = path+'reports/'
    os.makedirs(rep_path, exist_ok = True)
    out = open(rep_path+'response_qc_report.txt', 'w')

    # check #treatment in Monotherapy
    treat = list(set(monotherapy['.metadata_treatment']))
    print("Total treatment in monotherapy:", len(treat))
    out.write("Total treatment in monotherapy: "+str(len(treat))+'\n')
    #print(treat)

    # check #treatment in combination therapy
    treat = list(set(list(synergy['.metadata_treatment_1'])+list(synergy['.metadata_treatment_2'])))
    print("Total treatment in combination therapy:", len(treat))
    out.write("Total treatment in combination therapy: "+str(len(treat))+'\n')
    #print(treat)
    
    # check #cell line
    cl = list(set(monotherapy['.identifier_sample_name']))
    print("Total cell lines in monotherapy:", len(cl))
    out.write("Total cell lines in monotherapy: "+str(len(cl))+'\n')
    cl = list(set(synergy['.identifier_sample_name']))
    print("Total cell lines in combination therapy:", len(cl))
    out.write("Total cell lines in combination therapy: "+str(len(cl))+'\n')

    # check #tissue   
    tissue = list(set(monotherapy['.metadata_cancer_type']))
    print("Total cell lines in monotherapy:", len(tissue))
    out.write("Total cell lines in monotherapy: "+str(len(tissue))+'\n')
    tissue = list(set(synergy['.metadata_cancer_type']))
    print("Total cell lines in combination therapy:", len(tissue))
    out.write("Total cell lines in combination therapy: "+str(len(tissue))+'\n')

    # check #drug_pair tested in cell lines
    synergy['drug_pair'] = ['.'.join(sorted(list(r))) for _,r in synergy[['.metadata_treatment_1','.metadata_treatment_2']].iterrows()]
    cellline_summary = synergy.groupby('drug_pair')['.identifier_sample_name'].nunique().sort_values()
    mean = cellline_summary.mean()
    print("Averge # cellines tested for drug pairs:", mean)
    out.write("Averge # cellines tested for drug pairs: "+str(mean)+'\n')
    print(cellline_summary)
    cellline_summary.to_csv(rep_path+"drug_pair_cellline_summary.tsv", sep = '\t')

    # check #combinations
    synergy['combination'] = synergy['drug_pair']+'.'+synergy['.identifier_sample_name']
    combination_summary = synergy.groupby('combination')['combination'].count().sort_values()
    print(combination_summary)
    combination_summary.to_csv(rep_path+"combination_summary.tsv", sep = '\t')

    # get all drugs tested
    drugs_summary = monotherapy[['.metadata_treatment', '.metadata_moa']].drop_duplicates().sort_values(by = '.metadata_treatment').rename(columns = {'.metadata_treatment':"drug_name",  '.metadata_moa':"mode-of-action"})
    
    # prepare files for pulling drug target
    os.makedirs('../target_gene', exist_ok = True)
    if 'ho1' in fname:
        drugs_summary.to_csv("../target_gene/all_drugs_summary_ho1.csv", index = False)
    elif 'ho2' in fname:
        drugs_summary.to_csv("../target_gene/all_drugs_summary_ho2.csv", index = False)
    else:
        drugs_summary.to_csv("../target_gene/all_drugs_summary.csv", index = False) # will be used in drug target pull

    # check overall reprodicibility
    reproduce_path = rep_path+'reproducibility/'
    os.makedirs(reproduce_path, exist_ok = True)
    freq_sum, r1,p1,r2, p2 = check_reproducibility(responses, 'synergy', 'aoc')
    freq_sum.to_csv(reproduce_path+"response_all_replicates_aoc.tsv", sep ='\t', index = False )
    out.write("Overall reproducibility (AOC) Pearson: "+str(r1)+'\t'+str(p1)+ "; Spearman: "+str(r2)+'\t'+str(p2)+'\n')
    freq_sum, r1,p1,r2, p2 = check_reproducibility(responses, 'synergy', 'bliss')
    freq_sum.to_csv(reproduce_path+"response_all_replicates_bliss.tsv", sep ='\t', index = False )
    out.write("Overall reproducibility (Bliss) Pearson: "+str(r1)+'\t'+str(p1)+ "; Spearman: "+str(r2)+'\t'+str(p2)+'\n')

    # check monotherapy reprodicibility
    freq_sum, r1,p1,r2, p2 = check_reproducibility(monotherapy, 'monotherapy', 'aoc')
    freq_sum.to_csv(reproduce_path+"monotherapy_replicates_aoc.tsv", sep ='\t', index = False )
    out.write("Monotherapy reproducibility (AOC) Pearson: "+str(r1)+'\t'+str(p1)+ "; Spearman: "+str(r2)+'\t'+str(p2)+'\n')

    # check in-batch reprodicibility
    for batch in sorted(set(monotherapy['.identifier_batch'])):
        print(batch)
        monotherapy_new = monotherapy[monotherapy['.identifier_batch'] == batch]
        _,_,_,_,_ = check_reproducibility(monotherapy_new, 'monotherapy', 'aoc')

    # check synergy reprodicibility
    freq_sum, r1,p1,r2, p2 = check_reproducibility(synergy, 'synergy', 'aoc')
    freq_sum.to_csv(reproduce_path+"synergy_replictates_aoc.tsv", sep ='\t', index = False )
    out.write("Synergy reproducibility (AOC) Pearson: "+str(r1)+'\t'+str(p1)+ "; Spearman: "+str(r2)+'\t'+str(p2)+'\n')
    freq_sum, r1,p1,r2, p2 = check_reproducibility(synergy, 'synergy', 'bliss')
    freq_sum.to_csv(reproduce_path+"synergy_replicates_bliss.tsv", sep ='\t', index = False )
    out.write("Synergy reproducibility (Bliss) Pearson: "+str(r1)+'\t'+str(p1)+ "; Spearman: "+str(r2)+'\t'+str(p2)+'\n')
    
    # check in-batch reprodicibility
    for batch in sorted(set(synergy['.identifier_batch'])):
        print(batch)
        synergy_new = synergy[synergy['.identifier_batch'] == batch]
        _,_,_,_,_ = check_reproducibility(synergy_new, 'synergy', 'aoc')
        _,_,_,_,_ = check_reproducibility(synergy_new, 'synergy', 'bliss')
    
    out.close()

def chemical_structure_qc():
    """ QC of chemical structure information and feature preprocess

    Transform SMILE formatted chemical structure into 6 types of Fingerprints
    """
    # chemical structure required modules
    from rdkit import Chem,DataStructs
    from rdkit.Chem import MACCSkeys, AllChem
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit.Chem import rdmolops
    from pubchempy import get_compounds, get_cids, get_properties
    import openbabel
    import pybel

    def fp_to_feature(fp, max_bit):
        """
        """
        feature = [0 for i in range(max_bit)]
        for n in fp.bits:
            feature[n] = 1
        return feature

    df_map = pd.read_csv('../chemical_structure_features/merck_confidential_molecular_structures_20200420/smiles/smiles_and_ms_no.csv', header = 0)
    df_map['treatment_name_orig'] = df_map['treatment_name_orig'].replace('MSC[A-Z0-9]+', '', regex=True).replace('_+', '_', regex=True).replace('^_', '', regex=True).replace('_$', '', regex=True)
    df_map = df_map.drop_duplicates(subset = ['treatment_name_orig'])
    df_map = df_map.loc[df_map['treatment_name_orig'] != '']
    df_map= df_map.replace(updated_drug)
    print(df_map)

    # make new dirs
    os.makedirs("../chemical_structure_features/MACCS_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/Morgan_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/RDK_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/FP2_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/FP3_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/FP4_features/", exist_ok = True)

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
        np.savetxt('../chemical_structure_features/MACCS_features/'+r['treatment_name_orig'],np.array(feature_1))
        feature_1 = {'MACCS_'+str(idx+1):x for idx,x in enumerate(feature_1)}

        # Morgan Fingerprints (1024*1)
        fp = AllChem.GetMorganFingerprintAsBitVect(ms,2,nBits=1024)
        tmp = fp.ToBitString()
        feature_2 = list(map(int, tmp))
        np.savetxt('../chemical_structure_features/Morgan_features/'+r['treatment_name_orig'],np.array(feature_2))
        feature_2 = {'Morgan_'+str(idx+1):x for idx,x in enumerate(feature_2)}
        
        # FP2 (1024*1)
        fp = mol.calcfp('FP2')
        feature_3 = fp_to_feature(fp, max_bit_fp2)
        np.savetxt('../chemical_structure_features/FP2_features/'+r['treatment_name_orig'],np.array(feature_3))
        feature_3 = {'FP2_'+str(idx+1):x for idx,x in enumerate(feature_3)}

        # FP3 (56*1)
        fp = mol.calcfp('FP3')
        feature_4 = fp_to_feature(fp, max_bit_fp3)
        np.savetxt('../chemical_structure_features/FP3_features/'+r['treatment_name_orig'],np.array(feature_4))
        feature_4 = {'FP3_'+str(idx+1):x for idx,x in enumerate(feature_4)}

        # FP4 (308*1)
        fp = mol.calcfp('FP4')
        feature_5 = fp_to_feature(fp, max_bit_fp4)
        np.savetxt('../chemical_structure_features/FP4_features/'+r['treatment_name_orig'],np.array(feature_5))
        feature_5 = {'FP4_'+str(idx+1):x for idx,x in enumerate(feature_5)}

        # RDK Fingerprints (2048*1)
        fp = rdmolops.RDKFingerprint(ms)
        tmp = fp.ToBitString()
        feature_6 = list(map(int, tmp))
        np.savetxt('../chemical_structure_features/RDK_features/'+r['treatment_name_orig'],np.array(feature_6))
        feature_6 = {'RDK_'+str(idx+1):x for idx,x in enumerate(feature_6)}

        drug_name = {'treatment':r['treatment_name_orig']}
        the_features = {**drug_name,**feature_1,**feature_2,**feature_3,**feature_4,**feature_5,**feature_6}
        all_features.append(the_features)

    all_feature=pd.DataFrame.from_dict(all_features)
    all_feature.to_csv('../chemical_structure_features/all_chemical_structure.csv', index = False)


def cancer_dependency_qc():
    """ QC and feature construction of cancer dependency information
    # TODO: check this feature
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
    #response_qc(path = '../data/', fname = 'merck_confidential_responses_training_20200221.tsv')
    response_qc(path = '../data/hold_out/ho1/', fname = 'merck_confidential_responses_ho1_20200221.tsv')
    response_qc(path = '../data/hold_out/ho2/', fname = 'merck_confidential_responses_ho2_20200221.tsv')
    #chemical_structure_qc()
    #molecular_qc()
    #geneset_qc()
    #cancer_dependency_qc()

if __name__ == "__main__":
    main()

#!usr/bin/python3
#author:@rayezh
import numpy as np
import pandas as pd

data  = pd.read_csv('../merck_confidential_features_20200221.tsv', sep = '\t')
"""
Data types:

    1. Single gene molecular features: <GENESYMBOL>_snv/cnv/exp/lof;
    2. Geneset features: <GENESET>__coh_pat/lof_pat (geneset sep by '_');
    3. DDR features: *_ddr (78);

"""





################### split dataset ###############

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
print("uk:", len(uk), uk)
celltype = ['.metadata_cancer_type', '.metadata_cancer_subtype', '.identifier_sample_name']
##################### exp #######################

"""
#QC of gene expression level
"""

exp_cols = celltype+exp

exp_data = data[exp_cols]
exp_data.to_csv('exp_merck_confidential_features_20200221.tsv', sep = '\t', header = True)



##################### snv #######################
"""
#QC of number of alleles affected by pathigenic mutation or indel
"""

snv_cols = celltype+snv

snv_data = data[snv_cols]
snv_data.to_csv('snv_merck_confidential_features_20200221.tsv', sep = '\t', header = True)



##################### cnv #######################
"""
#QC of gene copy number variation
"""

cnv_cols = celltype+cnv

cnv_data = data[cnv_cols]
cnv_data.to_csv('cnv_merck_confidential_features_20200221.tsv', sep = '\t', header = True)


##################### lof #######################
"""
#QC of predicted gene loss of function
"""

lof_cols = celltype+lof

lof_data = data[cnv_cols]
lof_data.to_csv('lof_merck_confidential_features_20200221.tsv', sep = '\t', header = True)

#################### coh_pat ####################
coh_pat_cols = celltype+coh_pat

coh_pat_data = data[coh_pat_cols]
coh_pat_data.to_csv('coh_pat_merck_confidential_features_20200221.tsv', sep = '\t', header = True)
#################### lof_pat ####################
lof_pat_cols = celltype+lof_pat

lof_pat_data = data[lof_pat_cols]
lof_pat_data.to_csv('lof_pat_merck_confidential_features_20200221.tsv', sep = '\t', header = True)
##################### ddr #######################
ddr_cols = celltype + ddr

ddr_data = data[ddr_cols]
ddr_data.to_csv('ddr_merck_confidential_features_20200221.tsv', sep = '\t', header = True)








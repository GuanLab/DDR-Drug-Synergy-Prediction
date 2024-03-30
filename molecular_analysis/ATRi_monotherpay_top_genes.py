import pandas as pd
import numpy as np

all_exp = pd.read_csv('../feature/molecular_features/exp_merck_confidential_features_20200221.tsv', sep = '\t', index_col = 0)

mono = pd.read_csv('../feature/data/monotherapy_responses.tsv', sep = '\t')

#ATRi
mono_sub = mono.loc[mono['.metadata_moa']=="ATRi"]
print(mono_sub)

#all_genes = {'gene':[], '|cor|':[]}
all_genes = {'gene':[], 'cor':[]}
b = np.array(mono_sub['.response_aoc'])
for gene in all_exp.columns:
    if gene == 'ATR':
        continue
    a = np.array(all_exp[gene].loc[mono_sub['.identifier_sample_name']])
    cor = np.corrcoef(a,b)[0,1]
    all_genes['gene'].append(gene.split('_')[0])
    all_genes['cor'].append(cor)
    #all_genes['|cor|'].append(abs(cor))

#all_genes = pd.DataFrame.from_dict(all_genes).sort_values(by = '|cor|', ascending = False)
all_genes = pd.DataFrame.from_dict(all_genes).sort_values(by = 'cor', ascending = False)
print(all_genes)

all_genes.to_csv('ATRi_genes_all.tsv', sep = '\t', index = False)

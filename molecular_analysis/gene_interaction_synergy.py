import pandas as pd
import numpy as np
import json

'''
calculate bliss score between gene pairs based on expression level
'''

df = pd.read_csv('../feature/data/synergy_responses_with_monotherapy.tsv', sep = '\t', header = 0)
drug2gene = json.load(open('../feature/target_gene/drug2gene.json','r'))

interact_gene = {'gene1':[], 'gene2':[], 'efficacy':[], 'synergy':[]}

core_target = ['ATM', 'ATR', 'PRKDC']

for idx, row in df.iterrows():
    gene1s = drug2gene[row['.metadata_treatment_1']]
    gene2s = drug2gene[row['.metadata_treatment_2']]
    aoc = row['.response_aoc']
    bliss = row['.response_bliss']
    for g1 in gene1s:
        for g2 in gene2s:
            if g1 not in core_target:
                g3 = g1
                g1= g2
                g2 = g3
            interact_gene['gene1'].append(g1)
            interact_gene['gene2'].append(g2)
            interact_gene['efficacy'].append(aoc)
            interact_gene['synergy'].append(bliss)

interact_gene = pd.DataFrame.from_dict(interact_gene)
interact_gene =  interact_gene.groupby(['gene1', 'gene2'], as_index=False).mean()

df_aoc = []
df_bliss = []

for g in core_target:
    subdf_aoc = interact_gene.loc[interact_gene['gene1']==g].sort_values(by = ['efficacy'], ascending= False).iloc[:10,:3]
    df_aoc.append(subdf_aoc)
    subdf_bliss = interact_gene.loc[interact_gene['gene1']==g].sort_values(by = ['synergy'], ascending= False).iloc[:10,[0,1,3]]
    df_bliss.append(subdf_bliss)

df_aoc =pd.concat(df_aoc)
df_bliss = pd.concat(df_bliss)
df_aoc.to_csv('interaction_gene_aoc.tsv', sep = '\t', index = False, na_rep='NA')
df_bliss.to_csv('interaction_gene_bliss.tsv', sep = '\t', index = False, na_rep='NA')

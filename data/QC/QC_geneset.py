#!usr/bin/python3
import pandas as pd
import numpy as np
df = pd.read_csv('../merck_confidential_genesets_20200414.tsv', sep = '\t', header =0, index_col = 0)
print(len(list(set(df['geneset_name'].to_list())))) #9450
print(len(list(set(df['gene_symbol'].to_list())))) #18681

#read all target genes
allgenes = []
with open('../target_genes.txt', 'r') as f:
    for line in f:
        t = line.rstrip()
        allgenes.append(t)

allgenesets = list(set(df['geneset_name'].tolist()))
allgenesets = ','.join(allgenesets)
print(allgenesets)
f_geneset = open('all_geneset.txt','w')
f_geneset.write(allgenesets)
f_geneset.close()
#for every target gene, we make a 9550x1 vector, denoting the geneset feature for this gene
b_path = '../geneset_features/binary/'
w_path = '../geneset_features/weighted/'
#binray feature: 0 or 1
#weighted feature: 0~1 normalized all genes in this geneset
#make disctionary of weights(1/n, n is total genes in that geneset) of very geneset:
"""
w = {}
for i in allgenesets:
    r, c = df.loc[df['geneset_name'] == i].shape
    w.update({i:1/r})

for gene in allgenes:
    binary_feature = np.array([1 if x in df.loc[df['gene_symbol'] == gene]['geneset_name'].tolist() else 0 for x in allgenesets])
    np.savetxt(b_path+gene,binary_feature)
    weighted_feature = np.array([w[x] if x in df.loc[df['gene_symbol'] == gene]['geneset_name'].tolist() else 0 for x in allgenesets])
    np.savetxt(w_path+gene,weighted_feature)
"""

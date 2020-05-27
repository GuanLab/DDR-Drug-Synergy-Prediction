import pandas as pd
import pickle
import mygene
import math
from collections import defaultdict


"""
convert mouse gene-to-gene network to human
"""


human_genes = pickle.load(open('../all_genes.pkl','rb'))
h2m = pd.read_csv('HOM_MouseHumanSequence.rpt', sep = '\t', dtype = str, index_col = False)

#homologeneID


human = h2m.loc[h2m['Common Organism Name'] == 'human']


mouse = h2m.loc[h2m['Common Organism Name'] == 'mouse, laboratory']

entrez_2_mouse = defaultdict(lambda:[])
for idx, row in human.iterrows():
    homo = row['HomoloGene ID']
    homo_mouse = mouse.loc[mouse['HomoloGene ID'] == homo]
    entrez_2_mouse[row['EntrezGene ID']] = homo_mouse['Mouse MGI ID'].tolist()


mg = mygene.MyGeneInfo()
out = mg.querymany(human_genes, scopes='symbol', fields='entrezgene', species='human')

human_2_mouse = {}
for x in out:
    try:
        human_2_mouse.update({x['query']:entrez_2_mouse[x['entrezgene']]})
    except:
        pass

print(human_2_mouse)
print(len(human_2_mouse))
pickle.dump(human_2_mouse, open('human_to_mouse_dict.pkl','wb'))


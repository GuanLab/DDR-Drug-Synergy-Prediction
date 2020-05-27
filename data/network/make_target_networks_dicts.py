#!usr/bin/python
"""
make target gene networks and save as dictionaries
"""

import pandas as pd
import pickle


networks = pd.read_csv('network.tsv', sep ='\t', header = None)
target_genes = []
with open('../target_genes.txt', 'r') as f:
    for line in f:
        try:
            target = line.rstrip()
            targets = {}
            targets.update({row[1]:row[2] for idx,row in networks[networks[0]==target].iterrows()})
            targets.update({row[0]:row[2] for idx,row in networks[networks[1]==target].iterrows()})
            networkfile = './networks/'+target+'.pkl'
            pickle.dump(targets, open(networkfile, 'wb'))
        except:
            print(line)

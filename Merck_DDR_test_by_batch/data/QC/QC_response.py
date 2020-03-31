#!usr/bin/python3
#author: @rayezh
"""
QC of drug response data
"""
import pandas as pd
import numpy as np
import collections
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import scipy
from scipy import stats

def replicability_aoc(data):
    data = data[data['.response_aoc'].notna()]
    experiments = [frozenset(x) for _,x in data[['.identifier_sample_name', '.metadata_treatment_1','.metadata_treatment_2']].iterrows()]
    data['exp'] = experiments
    y = collections.Counter(experiments)
    rep_experiments = [{i:j} for i,j in y.items() if j > 1]
    rep = []
    freq = []
    for i in rep_experiments:
        rep.extend(list(i.keys()))
        freq.extend(list(i.values()))
    
    freq_summary = pd.DataFrame({"exp": rep, "freq": freq})  #replicates
    freq_summary.to_csv("syn_rep_freq_aoc.tsv", sep ='\t', index = False )
    
    cor_aoc =[]
    for rep in freq_summary['exp']:
        rows = data[data['exp'] == rep]
        r,c = rows.shape
        aoc = list(rows['.response_aoc'])
        cor_aoc.extend([[x,y] for x, y in combinations(aoc,2)])
        #print(cor_aoc)
        
    cor_aoc = np.array(cor_aoc)
    print(cor_aoc.shape)
    print(scipy.stats.spearmanr(cor_aoc[:,0], cor_aoc[:,1]))

def replicability_bliss(data):
    data = data[data['.response_bliss'].notna()]
    experiments = [frozenset(x) for _,x in data[['.identifier_sample_name', '.metadata_treatment_1','.metadata_treatment_2']].iterrows()]
    data['exp'] = experiments
    y = collections.Counter(experiments)
    rep_experiments = [{i:j} for i,j in y.items() if j > 1]
    rep = []
    freq = []
    for i in rep_experiments:
        rep.extend(list(i.keys()))
        freq.extend(list(i.values()))
    
    freq_summary = pd.DataFrame({"exp": rep, "freq": freq})  #replicates
    freq_summary.to_csv("syn_rep_freq.tsv", sep ='\t', index = False )
    
    cor_aoc =[]
    for rep in freq_summary['exp']:
        rows = data[data['exp'] == rep]
        r,c = rows.shape
        aoc = list(rows['.response_bliss'])
        cor_aoc.extend([[x,y] for x, y in combinations(aoc,2)])
        #print(cor_aoc)
        
        cor_aoc = np.array(cor_aoc)
        print(cor_aoc.shape)
        print(scipy.stats.pearsonr(cor_aoc[:,0], cor_aoc[:,1]))


data = pd.read_csv('../Synergy.tsv', sep = '\t')
#data = pd.read_csv('../Monotherapy.tsv', sep = '\t')
#replicability in Synvergy data
"""
replicability_bliss(data)
replicability_aoc(data)
"""












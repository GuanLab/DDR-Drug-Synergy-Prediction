#!usr/bin/python3
#author: @rayezh
import pandas as pd
import os
from sklearn.model_selection import train_test_split

features = pd.read_csv('merck_confidential_features_20200221.tsv', sep = '\t')
genes = pd.read_csv('merck_confidential_genesets_20200221.tsv', sep = '\t')
responses = pd.read_csv('merck_confidential_responses_training_20200221.tsv', sep = '\t')
responses = responses[responses['.metadata_treatment_1'].notna()]
responses = responses[responses['.metadata_treatment_2'].notna()]




#filter out all monotherapy
monotherapy = responses[responses['.metadata_treatment_1'] == "Monotherapy"].copy()
monotherapy = monotherapy.append(responses[responses['.metadata_treatment_2'] == "Monotherapy"])
#monotherapy = monotherapy.append(responses[responses['.metadata_treatment_2'].isna()]) for these data, aoc is also missing; might exclude(what is this?)
monotherapy.iloc[:,1:].to_csv('Monotherapy.tsv', sep ='\t', index = False,  na_rep = 'NA')

#construct feature (gene expression markers) matrix




#Synergy files(for 5-fold cross validation)
synergy = responses[responses['.metadata_treatment_1'] != "Monotherapy"].copy()
synergy = synergy[synergy['.metadata_treatment_2'] != "Monotherapy"].copy()
synergy = synergy[synergy['.metadata_treatment_2'].notna()].copy()
synergy.iloc[:,1:].to_csv('Synergy.tsv', sep ='\t', index = False,  na_rep = 'NA')
synergy = synergy.iloc[:,1:]


def test_train_split_by_cell_line(synergy):
    """
    #test set 1: test within indications
    *split by cell lines(balanced by indications)
    """
    cell_lines = list(set(list(synergy['.identifier_sample_name'])))
    print(sorted(cell_lines))
    print(len(cell_lines))
    
    #5-fold cross validation
    for seed in range(5):
        filepath = '../test_by_cell_line/fold_'+str(seed)
        os.mkdir(filepath)
        train,test = train_test_split(cell_lines, shuffle = True, test_size = 0.2,random_state = seed)
        train_df = synergy[synergy['.identifier_sample_name'].isin(train)]
        test_df = synergy[synergy['.identifier_sample_name'].isin(test)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)

def test_train_split_by_indications(synergy):
    """
    test set 2:test in a new indication
    *split by indication (Indication means tissure? cancer type? subtype?)
    
    Question:
    1. the Merck hold-out set 2 contains one single indication (of all cell lines);
    So, how should we test for the cross cell line performance?
    My idea is: 
    1. we cross test all indications by leave-one-out and see average.
    """

    indications = [frozenset((x['.metadata_treatment_1'], x['.metadata_treatment_2'])) for _,x in synergy.iterrows()]
    indications = list(set(indications))
    print(indications)
    print(len(indications))

    for seed in range(5):
        filepath = '../test_by_indication/fold_'+str(seed)
        os.mkdir(filepath)
        train,test = train_test_split(indications, shuffle = True, random_state = seed)
        train_df = synergy[frozenset(synergy['.metadata_treatment_1'],synergy['.metadata_treatment_2']).isin(train)]
        test_df = synergy[synergy['.metadata_treatment_1'].isin(test) & synergy['.metadata_treatment_2'].isin(test)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)


def test_train_split_by_batch(synergy):
    """
    test set *: independent batch
    *split by batch
    """

    batches = list(set(list(synergy['.identifier_batch'])))
    print(sorted(batches))
    print(len(batches))
    
    for seed in range(5):
        filepath = '../test_by_batch/fold_'+str(seed)
        os.mkdir(filepath)
        train,test = train_test_split(batches, shuffle = True, test_size = 0.25,random_state = seed)
        train_df = synergy[synergy['.identifier_batch'].isin(train)]
        test_df = synergy[synergy['.identifier_batch'].isin(test)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)




test_train_split_by_batch(synergy)








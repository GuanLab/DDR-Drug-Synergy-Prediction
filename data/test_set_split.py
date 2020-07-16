#!usr/bin/python3
#author: @rayezh
import pandas as pd
import os
from sklearn.model_selection import train_test_split, KFold

features = pd.read_csv('merck_confidential_features_20200221.tsv', sep = '\t')
#genes = pd.read_csv('merck_confidential_genesets_20200221.tsv', sep = '\t')
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
    os.makedirs('../test_by_cell_line/', exist_ok = True)
    kf = KFold(n_splits = 5)
    i = 0
    #5-fold cross validation
    for train, test in kf.split(cell_lines):
        filepath = '../test_by_cell_line/fold_'+str(i)
        i+=1
        os.makedirs(filepath, exist_ok = True)
        train_f = [cell_lines[j] for j in train]
        test_f = [cell_lines[j] for j in test]
        train_df = synergy[synergy['.identifier_sample_name'].isin(train_f)]
        test_df = synergy[synergy['.identifier_sample_name'].isin(test_f)]
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
    os.makedirs('../test_by_batch/', exist_ok = True)
    for seed in range(5):
        filepath = '../test_by_batch/fold_'+str(seed)
        os.makedirs(filepath, exist_ok = True)
        train,test = train_test_split(batches, shuffle = True, test_size = 0.25,random_state = seed)
        train_df = synergy[synergy['.identifier_batch'].isin(train)]
        test_df = synergy[synergy['.identifier_batch'].isin(test)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)

def test_train_split_by_indication(synergy):
    """ 
    test set *: independent indication (tissue)
    *split by indication(tissue)
    train on one tissue and test on the remaining 11 tissues
    """
    cell_lines = list(set(list(synergy['.metadata_cancer_type'])))
    os.makedirs('../test_by_indication/', exist_ok = True)
    #tissue specific train/validation
    for train in cell_lines:
        train_f = [train]
        test_f = cell_lines.copy()
        test_f.remove(train)
        filepath = '../test_by_indication/fold_'+train
        os.makedirs(filepath, exist_ok = True)
        #print(train_f)
        #print(test_f)
        train_df = synergy[synergy['.metadata_cancer_type'].isin(train_f)]
        test_df = synergy[synergy['.metadata_cancer_type'].isin(test_f)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)

def split_moa(f):
    import glob
    fpaths = glob.glob(f+'*')
    for path in fpaths:
        os.makedirs(path+'/moa', exist_ok = True)
        df = pd.read_csv(path+'/Test.tsv', sep = '\t')
        moa_pairs = ['_'.join(sorted([r['.metadata_moa_1'],r['.metadata_moa_2']])) for _,r in df.iterrows()]
        df['moa_pairs'] = moa_pairs
        moa_pairs = list(set(moa_pairs))
        for p in moa_pairs:
            df_new = df.loc[df['moa_pairs'] == p]
            df_new.to_csv(f+'/moa/'+p+'.tsv', index = False, sep = '\t')



test_train_split_by_cell_line(synergy)
test_train_split_by_batch(synergy)
test_train_split_by_indication(synergy)

split_moa('../test_by_cell_line/')
split_moa('../test_by_batch/')
split_moa('../test_by_indication/')






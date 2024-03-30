#!usr/bin/python3
#author: @rayezh
import pandas as pd
import json,os
from sklearn.model_selection import train_test_split, KFold, StratifiedKFold

def test_train_split_by_cell_line(synergy):
    """
    #test set 1: test within indications
    *split by cell lines(balanced by indications)
    """
    tmp = synergy[['.identifier_sample_name', '.metadata_cancer_type']].drop_duplicates().sort_values(by = '.identifier_sample_name').reset_index(drop = True)
    #print(tmp)
    cell_lines = list(tmp['.identifier_sample_name'])
    cancer_types = list(tmp['.metadata_cancer_type'])
    print(cell_lines, len(cell_lines))
    print(cancer_types, len(cancer_types))
    os.makedirs('../../test_by_cell_line/', exist_ok = True)
    # stratefied fold: each fold has the same proportion of cell lines from each indication
    # kf = KFold(n_splits = 5, shuffle = True, random_state = 42)
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=0)
    i = 0
    #5-fold cross validation
    for train, test in skf.split(cell_lines, cancer_types):
        filepath = '../../test_by_cell_line/fold_'+str(i)
        i+=1
        os.makedirs(filepath, exist_ok = True)
        train_f = [cell_lines[j] for j in train]
        test_f = [cell_lines[j] for j in test]
        #print(test_f, len(set([cancer_types[j] for j in test])), len(set([cancer_types[j] for j in train])))
        train_df = synergy[synergy['.identifier_sample_name'].isin(train_f)]
        test_df = synergy[synergy['.identifier_sample_name'].isin(test_f)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)
    
    # put hold-out set ho1 in the validation set
    fpath = '../../test_by_cell_line/hold_out_validation/'
    os.makedirs(fpath, exist_ok = True)
    val_data = pd.read_csv('../data/hold_out/ho1/synergy_responses_with_monotherapy.tsv', sep ='\t')
    val_data.to_csv(fpath+'/Test.tsv', sep = '\t', index = False)

def test_train_split_by_batch(synergy):
    """
    test set *: independent batch
    *split by batch
    """

    batches = sorted(list(set(list(synergy['.identifier_batch']))))
    main_batch = ["Oncolead_007", "Oncolead_010", "Oncolead_012"]
    print(batches)
    print(len(batches))
    os.makedirs('../../test_by_batch/', exist_ok = True)
    i = 0
    seed  = 0
    while i < 5:
        seed += 1
        train_f,test_f = train_test_split(batches, shuffle = True, test_size = 0.25,random_state = seed)
        if len(list(set(train_f) & set(main_batch))) == 2:
            print(train_f)
            print(test_f)
            filepath = '../../test_by_batch/fold_'+str(i)
            os.makedirs(filepath, exist_ok = True)
            train_df = synergy[synergy['.identifier_batch'].isin(train_f)]
            test_df = synergy[synergy['.identifier_batch'].isin(test_f)]
            train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
            test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)
            i+=1

def test_train_split_by_indication(synergy):
    """
    leave-one-out cross validation 
    test set *: independent cancer type
    *split by cancer type
    """
    indications = sorted(list(set(list(synergy['.metadata_cancer_type']))))
    print(indications)
    print(len(indications))
    os.makedirs('../../test_by_indication/', exist_ok = True)
    i = 0
    for test_f in indications:
        train_f = indications.copy()
        train_f.remove(test_f)
        filepath = '../../test_by_indication/fold_'+test_f
        os.makedirs(filepath, exist_ok = True)
        i += 1
        print(test_f)
        train_df = synergy[synergy['.metadata_cancer_type'].isin(train_f)]
        test_df = synergy[synergy['.metadata_cancer_type'].isin([test_f])]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)
    
    # put hold-out set ho2 in the validation set
    fpath = '../../test_by_indication/hold_out_validation/'
    os.makedirs(fpath, exist_ok = True)
    val_data = pd.read_csv('../data/hold_out/ho2/synergy_responses_with_monotherapy.tsv', sep ='\t')
    indications = sorted(list(set(list(val_data['.metadata_cancer_type']))))
    for test_i in indications:
        test = val_data[val_data['.metadata_cancer_type'] == test_i]
        test.to_csv(fpath+'/Test_'+test_i+'.tsv', sep = '\t', index = False)



def test_train_split_by_moa_pair(synergy):
    """ 
    test set *: independent moa pairs
    *split by moa pair
    """
    # TODO: new moa-known moa; new moa-new moa; known moa-known moa
    moa_pairs = ['_'.join(sorted([r['.metadata_moa_1'],r['.metadata_moa_2']])) for _,r in synergy.iterrows()]
    synergy['moa_pair'] = moa_pairs
    print(moa_pairs)
    print(len(moa_pairs))
    os.makedirs('../../test_by_moa_pairs/', exist_ok = True)
    kf = KFold(n_splits = 5, shuffle = True, random_state = 42) 
    i = 0 
    for train, test in kf.split(moa_pairs):
        filepath = '../../test_by_moa_pairs/fold_'+str(i)
        os.makedirs(filepath, exist_ok = True)
        i += 1
        train_f = [moa_pairs[j] for j in train]
        test_f = [moa_pairs[j] for j in test]
        print(test_f)
        train_df = synergy[synergy['moa_pair'].isin(train_f)]
        test_df = synergy[synergy['moa_pair'].isin(test_f)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)


def test_train_split_by_drug_pair(synergy):
    """
    test set *: independent drug pairs
    *split by drug pair
    """
    # TODO: new drug-known drug; new drug-new drug; known drug-known drug
    moa_pairs = ['_'.join(sorted([r['.metadata_treatment_1'],r['.metadata_treatment_2']])) for _,r in synergy.iterrows()]
    synergy['treatment_pair'] = moa_pairs
    print(moa_pairs)
    print(len(moa_pairs))
    os.makedirs('../../test_by_treatment_pairs/', exist_ok = True)
    kf = KFold(n_splits = 5, shuffle = True, random_state = 42)
    i = 0
    for train, test in kf.split(moa_pairs):
        filepath = '../../test_by_treatment_pairs/fold_'+str(i)
        os.makedirs(filepath, exist_ok = True)
        i += 1
        train_f = [moa_pairs[j] for j in train]
        test_f = [moa_pairs[j] for j in test]
        print(test_f)
        train_df = synergy[synergy['treatment_pair'].isin(train_f)]
        test_df = synergy[synergy['treatment_pair'].isin(test_f)]
        train_df.to_csv(filepath+'/Train.tsv', sep = '\t', index = False)
        test_df.to_csv(filepath+'/Test.tsv', sep = '\t', index = False)

def test_train_split_by_indication_specific(synergy):
    """ 
    test set *: independent indication (tissue)
    *split by indication(tissue)
    train on one tissue and test on the remaining 11 tissues
    """
    cancer_type = sorted(list(set(list(synergy['.metadata_cancer_type']))))
    os.makedirs('../../test_by_indication_specific/', exist_ok = True)
    #tissue specific train/validation
    for train in cancer_type:
        train_f = [train]
        test_f = cancer_type.copy()
        test_f.remove(train)
        filepath = '../../test_by_indication_specific/fold_'+train
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
        print(path)
        os.makedirs(path+'/moa', exist_ok = True)
        df = pd.read_csv(path+'/Test.tsv', sep = '\t')
        moas = ["ATMi", "ATRi", "DNAPKi"]
        moa_pairs = ['_'.join(sorted([r['.metadata_moa_1'],r['.metadata_moa_2']])) for _,r in df.iterrows()]
        #df['moa_pairs'] = moa_pairs
        #moa_pairs = list(set(moa_pairs))
        for m in moas:
            df_new = df.loc[(df[['.metadata_moa_1','.metadata_moa_2']] == m).any(axis = 1)]
            df_new.to_csv(path+'/moa/'+m+'.tsv', index = False, sep = '\t')

def split_all_moa(f):
    import glob
    fpaths = glob.glob(f+'*')
    for path in fpaths:
        print(path)
        os.makedirs(path+'/all_moa', exist_ok = True)
        df = pd.read_csv(path+'/Test.tsv', sep = '\t')
        moa_pairs = ['_'.join(sorted([r['.metadata_moa_1'],r['.metadata_moa_2']])) for _,r in df.iterrows()]
        df['moa_pairs'] = moa_pairs
        moa_pairs = list(set(moa_pairs))
        for p in moa_pairs:
            df_new = df.loc[df['moa_pairs'] == p]
            df_new.to_csv(path+'/all_moa/'+p+'.tsv', index = False, sep = '\t')

def split_cell_line(f):
    import glob
    fpaths = glob.glob(f+'*')
    for path in fpaths:
        os.makedirs(path+'/cell_line', exist_ok = True)
        df = pd.read_csv(path+'/Test.tsv', sep = '\t')
        cell_lines = set(list(df['.identifier_sample_name']))
        for c in cell_lines:
            df_new = df.loc[df['.identifier_sample_name'] == c]
            df_new.to_csv(path+'/cell_line/'+c+'.tsv', index = False, sep = '\t')


def split_tissue(f):
    import glob
    fpaths = glob.glob(f+'*')
    for path in fpaths:
        os.makedirs(path+'/tissue', exist_ok = True)
        df = pd.read_csv(path+'/Test.tsv', sep = '\t')
        tissue = set(list(df['.metadata_cancer_type']))
        for c in tissue:
            df_new = df.loc[df['.metadata_cancer_type'] == c]
            df_new.to_csv(path+'/tissue/'+c+'.tsv', index = False, sep = '\t')


syn_data = pd.read_csv('../data/synergy_responses_with_monotherapy.tsv', sep ='\t') # use synergy without missing monotherapy 

test_train_split_by_cell_line(syn_data)
test_train_split_by_indication(syn_data)
#test_train_split_by_moa_pair(syn_data)
#test_train_split_by_drug_pair(syn_data)
#test_train_split_by_batch(syn_data)

#split_moa('../../test_by_cell_line/')
#split_moa('../../test_by_indication/')

#split_all_moa('../../test_by_cell_line/')
#split_all_moa('../../test_by_indication/')

#split_tissue('../../test_by_indication/')
#split_cell_line('../../test_by_cell_line/')

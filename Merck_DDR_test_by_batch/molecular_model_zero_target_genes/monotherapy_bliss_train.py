#!usr/bin/python3
#author:@rayezh
import pandas as pd
import numpy as np
import pickle
import json

def train_RandomForest_model(X, Y):
    from sklearn.ensemble import RandomForestRegressor
    regressor = RandomForestRegressor(n_estimators=500, random_state=0)
    regressor.fit(X, Y)
    return regressor

def train_LightGBM_model(X,Y):
    import lightgbm as lgb
    param = {'boosting_type': 'gbdt',
            'objective': 'regression',
            'num_leaves': 10,
            'learning_rate': 0.05,
            'verbose': 0,
            'n_estimators': 300,
            'reg_alpha': 2.0,
                   }
    train_data = lgb.Dataset(np.asarray(X), np.asarray(Y),categorical_feature=[0,1,2,3,30,31,32,33])
    regressor =lgb.train(param, train_data, num_boost_round=1000)
    return regressor

def pd_to_XY(data, if_train = True):
    target_genes = []
    #read moa's target genes
    f = open('../data/moa_to_genes.json','r').read()
    moa2gene = json.loads(f)
    
    gene_features = pd.read_csv('../data/merck_confidential_features_20200221.tsv', sep ='\t', index_col = 0)
    cols = gene_features.columns.to_list()

    #quantile:
    batches = list(set(list(data['.identifier_batch'])))
    mean_ba = {batch:np.mean(np.array(data.loc[data['.identifier_batch'] == batch]['.response_bliss'])) for batch in batches}
    std_ba = {batch:np.std(np.array(data.loc[data['.identifier_batch'] == batch]['.response_bliss'])) for batch in batches}
    
    embed = {}
    for i in range(13):
        embed.update({'Oncolead_'+'{:0>3}'.format(i):np.array([i])})
    X = []
    Y = []
    for _, row in data.iterrows():
        try:
            #feature0 =embed[row['.identifier_batch']]
            feature1 = np.loadtxt('./cell_drug_features/'+row['.identifier_sample_name']+'.'+row['.metadata_treatment_1'])
            feature2 = np.loadtxt('./cell_drug_features/'+row['.identifier_sample_name']+'.'+row['.metadata_treatment_2'])
            
            target_genes = list(set(moa2gene[row['.metadata_moa_1']]+moa2gene[row['.metadata_moa_2']])) #overlapped target genes of two drugs
            cols = gene_features.columns.to_list()
            tmp = [list(filter(lambda x: subs in x, cols)) for subs in target_genes]
            target_genes = []
            for x in tmp:
                target_genes.extend(x)
            #print(len(target_genes))
            mol = gene_features.loc[gene_features['.identifier_sample_name'] == row['.identifier_sample_name']].copy()
            mol[target_genes] = 0
            feature3 = np.array(mol.iloc[0,3:].to_list())
            #print(feature3)
            x = np.concatenate((feature1,feature2, feature3))
            #print(x.shape)
            y = row['.response_bliss']
            X.append(x)
            if if_train:
                Y.append((y-mean_ba[row['.identifier_batch']])/std_ba[row['.identifier_batch']])
            else:
                Y.append(y)
        except:
            pass
    
    return X, Y

############# five fold cross-validation ###################

results = open('cor_bliss.txt', 'w')
for i in range(5):
    train_path = '../test_by_batch/fold_'+str(i)+'/Train.tsv'
    test_path = '../test_by_batch/fold_'+str(i)+'/Test.tsv'

    #construct synbergy features from two monotherapy features
    Train = pd.read_csv(train_path, sep = '\t')
    #predict aoc, exclude those with nan in aoc
    Train = Train[Train['.response_bliss'].notna()]
    Train_X,Train_Y = pd_to_XY(Train,  if_train = True)
    print(len(Train_X), len(Train_Y))
    regressor = train_LightGBM_model(Train_X, Train_Y)
    print('saving model ...')
    filename = 'monotherapy_bliss'+str(i)+'.model'
    pickle.dump(regressor, open(filename, 'wb'))
    
    beta_Y = regressor.predict(Train_X)
    cor_beta = np.corrcoef(beta_Y, Train_Y)[0,1]
    print(cor_beta)

    Test = pd.read_csv(test_path, sep = '\t')
    #predict aoc, exclude those with nan in aoc
    Test = Test[Test['.response_bliss'].notna()]
    Test_X,Test_Y = pd_to_XY(Test, if_train = False)
    print(len(Test_X), len(Test_Y))
    print('start prediction ...')
    pred_Y = regressor.predict(Test_X)

    cor  = np.corrcoef(pred_Y, Test_Y)[0,1]
    print('correlation =',cor)
    results.write('%.4f\n'%(cor))

results.close()

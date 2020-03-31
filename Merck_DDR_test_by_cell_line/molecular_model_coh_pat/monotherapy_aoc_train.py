#!usr/bin/python3
#author:@rayezh
import pandas as pd
import numpy as np
import pickle


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
            'n_estimators': 800,
            'reg_alpha': 2.0,
                   }
    train_data = lgb.Dataset(np.asarray(X), np.asarray(Y),categorical_feature=[0,1,2,3,4,31,32,33,34])
    regressor =lgb.train(param, train_data, num_boost_round=1000)
    return regressor

def pd_to_XY(data):
    gene_features = pd.read_csv('../data/coh_pat_merck_confidential_features_20200221.tsv', sep ='\t', index_col = 0)
    embed = {}
    for i in range(13):
        embed.update({'Oncolead_'+'{:0>3}'.format(i):np.array([i])})
    X = []
    Y = []
    for _, row in data.iterrows():
        try:
            feature0 =embed[row['.identifier_batch']]
            feature1 = np.loadtxt('./cell_drug_features/'+row['.identifier_sample_name']+'.'+row['.metadata_treatment_1'])
            feature2 = np.loadtxt('./cell_drug_features/'+row['.identifier_sample_name']+'.'+row['.metadata_treatment_2'])
            feature3 = np.array(gene_features.loc[gene_features['.identifier_sample_name'] == row['.identifier_sample_name']].iloc[0,3:].to_list())
            x = np.concatenate((feature0,feature1,feature2, feature3))
            y = row['.response_aoc']
            X.append(x)
            Y.append(y)
        except:
            pass
    
    return X, Y

############# five fold cross-validation ###################

results = open('cor_aoc.txt', 'w')
for i in range(5):
    train_path = '../test_by_cell_line/fold_'+str(i)+'/Train.tsv'
    test_path = '../test_by_cell_line/fold_'+str(i)+'/Test.tsv'

    #construct synbergy features from two monotherapy features
    Train = pd.read_csv(train_path, sep = '\t')
    #predict aoc, exclude those with nan in aoc
    Train = Train[Train['.response_aoc'].notna()]
    Train_X,Train_Y = pd_to_XY(Train)
    print(len(Train_X), len(Train_Y))
    regressor = train_LightGBM_model(Train_X, Train_Y)
    print('saving model ...')
    filename = 'monotherapy_aoc'+str(i)+'.model'
    pickle.dump(regressor, open(filename, 'wb'))
    
    beta_Y = regressor.predict(Train_X)
    cor_beta = np.corrcoef(beta_Y, Train_Y)[0,1]
    print(cor_beta)

    Test = pd.read_csv(test_path, sep = '\t')
    #predict aoc, exclude those with nan in aoc
    Test = Test[Test['.response_aoc'].notna()]
    Test_X,Test_Y = pd_to_XY(Test)
    print(len(Test_X), len(Test_Y))
    print('start prediction ...')
    pred_Y = regressor.predict(Test_X)

    cor  = np.corrcoef(pred_Y, Test_Y)[0,1]
    print('correlation =',cor)
    results.write('%.4f\n'%(cor))

results.close()

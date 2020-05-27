#!usr/bin/python3
#author:@rayezh
import pandas as pd
import numpy as np
import pickle
import json
from tqdm import tqdm
import lightgbm as lgb
import shap
import matplotlib.pyplot as plt
import glob
import os
from utils import pd_to_XY

def train_RandomForest_model(X, Y):
    from sklearn.ensemble import RandomForestRegressor
    regressor = RandomForestRegressor(n_estimators=500, random_state=0)
    regressor.fit(X, Y)
    return regressor

def train_LightGBM_model(X,Y):
    import lightgbm as lgb
    param = {'boosting_type': 'gbdt',
            'objective': 'regression',
            'num_leaves': 20,
            'learning_rate': 0.05,
            'verbose': 0,
            'n_estimators': 1000,
            'reg_alpha': 2.0,
                   }
    train_data = lgb.Dataset(X, Y,categorical_feature=[0,1,2,3,4,5,4659,4660,4661,4662,4663])
    regressor =lgb.train(param, train_data, num_boost_round=1000)
    return regressor



def train_models(target):
    results = open('cor_'+target+'.txt', 'w')
    for i in range(5):
        train_path = '../test_by_cell_line/fold_'+str(i)+'/Train.tsv'
        test_path = '../test_by_cell_line/fold_'+str(i)+'/Test.tsv'
        
        #construct synbergy features from two monotherapy features
        Train = pd.read_csv(train_path, sep = '\t')
        #predict aoc, exclude those with nan in aoc
        Train = Train[Train['.response_'+target].notna()]
        Train_X,Train_Y =pd_to_XY(target, Train, if_train = True)
        
        print('start training model ...')
        regressor = train_LightGBM_model(Train_X, Train_Y)
        print('saving model ...')
        filename = 'monotherapy_'+target+'_'+str(i)+'.model'
        pickle.dump(regressor, open(filename, 'wb'))
        
        beta_Y = regressor.predict(Train_X)
        cor_beta = np.corrcoef(beta_Y, Train_Y)[0,1]
        print(cor_beta)
        
        Test = pd.read_csv(test_path, sep = '\t')
        #predict aoc, exclude those with nan in aoc
        Test = Test[Test['.response_'+target].notna()]
        Test_X,Test_Y = pd_to_XY(target, Test, if_train = False)
        #print(len(Test_X), len(Test_Y))
        print('start prediction ...')
        pred_Y = regressor.predict(Test_X)
        
        cor  = np.corrcoef(pred_Y, Test_Y)[0,1]
        
        print('correlation =',cor)
        results.write('%.4f\n'%(cor))
        shap_values = shap.TreeExplainer(regressor).shap_values(Test_X)
        shap.summary_plot(shap_values, Test_X, show=False)
        plt.savefig('importance_scatter_'+target+'_'+str(i)+'_xgb.pdf', format='pdf', dpi=1200, bbox_inches='tight')
        plt.close()
        
        shap_df =pd.DataFrame(np.array(shap_values), columns = Test_X.columns.tolist())
        shap_df.to_csv('SHAP_'+target+'_'+str(i)+'.csv', index = False)
    results.close()

def SHAP_analysis_moa(target):
    for i in range(5):
        filename = 'monotherapy_'+target+'_'+str(i)+'.model'
        regressor = pickle.load(open(filename, 'rb'))
        test_path = '../test_by_cell_line/fold_'+str(i)+'/moa/*.tsv'
        all_moa_path = glob.glob(test_path)
        os.makedirs('moa', exist_ok = True)
        for path in all_moa_path:
            try:
                name = path.split('/')[-1].split('.')[0]
                os.makedirs('./moa/'+name, exist_ok = True)
                results = open('./moa/'+name+'/cor_'+target+'.txt', 'a')
                Test = pd.read_csv(path, sep = '\t')
                #predict aoc, exclude those with nan in aoc
                Test = Test[Test['.response_'+target].notna()]
                Test_X,Test_Y = pd_to_XY(target, Test, if_train = False)
                #print(len(Test_X), len(Test_Y))
                print('start prediction ...')
                pred_Y = regressor.predict(Test_X)
        
                cor  = np.corrcoef(pred_Y, Test_Y)[0,1]
                print('correlation =',cor)
                results.write('%.4f\n'%(cor))
                shap_values = shap.TreeExplainer(regressor).shap_values(Test_X)
                shap.summary_plot(shap_values, Test_X, show=False)
                plt.savefig('./moa/'+name+'/importance_scatter_'+target+'_'+str(i)+'_xgb.pdf', format='pdf', dpi=1200, bbox_inches='tight')
                plt.close()
        
                shap_df =pd.DataFrame(np.array(shap_values), columns = Test_X.columns.tolist())
                shap_df.to_csv('./moa/'+name+'/SHAP_'+target+'_'+str(i)+'.csv', index = False)
                results.close()
            except:
                print(path)

def Dependency_plot_of_genes(target,moa ='ATMi_ATRi', gene = 'PARP1'):
    for i in range(5):
        try:

            filename = 'monotherapy_bliss_'+str(i)+'.model'
            regressor = pickle.load(open(filename, 'rb'))
            moa_path = '../test_by_cell_line/fold_'+str(i)+'/'+moa+'.tsv'

            Test = pd.read_csv(moa_path, sep = '\t')
            Test = Test[Test['.response_'+target].notna()]
            Test_X,Test_Y = pd_to_XY(target, Test, if_train = False)

            print('start prediction ...')
            pred_Y = regressor.predict(Test_X)

            cor  = np.corrcoef(pred_Y, Test_Y)[0,1]
            print('correlation =',cor)
            shap_values = shap.TreeExplainer(regressor).shap_values(Test_X)
            shap.dependence_plot(gene+"_exp", shap_values, Test_X, interaction_index=gene+"_cnv")
            plt.savefig('dependency_'+gene+'_exp_cnv_'+target+'_'+str(i)+'_lgb.pdf', format='pdf', dpi=1200, bbox_inches='tight')
            plt.close()

            shap.dependence_plot(gene+"_exp", shap_values, Test_X, interaction_index=gene+"_snv")
            plt.savefig('dependency_'+gene+'_exp_snv_'+target+'_'+str(i)+'_lgb.pdf', format='pdf', dpi=1200, bbox_inches='tight')
            plt.close()

            shap.dependence_plot(gene+"_exp", shap_values, Test_X, interaction_index=gene+"_lof")
            plt.savefig('dependency_'+gene+'_exp_lof_'+target+'_'+str(i)+'_lgb.pdf', format='pdf', dpi=1200, bbox_inches='tight')
            plt.close()

        except:
            print('fail to make dependency plot.')





train_models('aoc')
train_models('bliss')

SHAP_analysis_moa('aoc')
SHAP_analysis_moa('bliss')

Dependency_plot_of_genes('aoc',moa ='ATMi_ATRi', gene = 'PARP1')
Dependency_plot_of_genes('bliss',moa ='ATMi_ATRi', gene = 'PARP1')

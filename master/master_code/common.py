#!usr/bin/python3
#author:@rayezh
import pandas as pd
import numpy as np
import numpy.ma as ma
from scipy import sparse 
from tqdm import tqdm
import lightgbm as lgb
import shap
import random
import matplotlib.pyplot as plt
from glob import glob
import pickle, json, os
from build_feature_dataset import *
from utils import *
from shap_analysis import *
from models import *

def train(train_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path, surrogate_gene, surrogate_geneset, surrogate_chem):
    """ Train LightGBM synergy prediction model
    
    Params
    ------
    train_path: str
        path to training set
    target: str
        prediction target(aoc or bliss)
    model_path: str
        path where to save the trained  model

    Yields
    ------
    cor: float
        correlation between prediction and gold standard on training dataset
    """
    print("Loading Training dataset ...")
    Train = pd.read_csv(train_path, sep = '\t')
        
    # Exclude datasets with missing target score
    Train = Train[Train['.response_'+target].notna()]

    # split validation set by half cell lines/tissue
    #random.seed(0)
    #if 'test_by_cell_line' in train_path:
    #    cell_lines = sorted(set(Train['.identifier_sample_name']))
    #    cell_lines_val = random.sample(cell_lines,int(0.5*len(cell_lines)))
    #    Val = Train.loc[Train['.identifier_sample_name'].isin(cell_lines_val),:]
    #    Train = Train.loc[~Train['.identifier_sample_name'].isin(cell_lines_val),:]
    #elif 'test_by_cell_line' in train_path:
    #    tissues = sorted(set(Train['.metadata_cancer_type']))
    #    tissues_val = random.sample(tissues,int(0.5*len(tissues)))
    #    Val = Train.loc[Train['.metadata_cancer_type'].isin(tissues_val),:]
    #    Train = Train.loc[~Train['.metadata_cancer_type'].isin(tissues_val),:]
        
    # Build Training feature set
    f_name, Train_X,Train_Y = build_feature_dataset(Train, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem, if_train = True)
    #f_name, Val_X, Val_Y = build_feature_dataset(Val, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem, if_train = True)    
    
    print('start training model ...')
    #regressor = train_LightGBM_model(Train_X, Train_Y,  Val_X, Val_Y, f_name)
    regressor = train_LightGBM_model(Train_X, Train_Y, f_name)
    print('saving model ...')
    pickle.dump(regressor, open(model_path, 'wb'))
        
    #beta_Y = regressor.predict(Train_X, num_iteration=regressor.best_iteration)
    beta_Y = regressor.predict(Train_X)
    cor = np.corrcoef(beta_Y, Train_Y)[0,1]
    print("Prediction-gold standard's correlation on training set: ", cor)
    
    return cor

def predict(test_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path, surrogate_gene, surrogate_geneset, surrogate_chem):
    """ Perform trained models on test dataset
    
    Params
    ------
    test_path: str
        path to test set
    target: str
        prediction target(aoc or bliss)
    model_path: str
        path to saved model

    Yields
    ------
    cor: float
        correlation between prediction and gold standard on training dataset
    pred: Numpy array
        1st column: gold standard
        2nd column: prediction
    """
    print("Loading Test dataset ...")
    Test = pd.read_csv(test_path, sep = '\t')

    # Exclude datasets with missing target score
    #Test = Test[Test['.response_'+target].notna()]
    if Test.shape[0] > 0:
        f_name, Test_X, Test_Y = build_feature_dataset(Test, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem, if_train = False)
        print('Loading saved model ...')
        regressor = pickle.load(open(model_path, 'rb'))
        
        print('start prediction ...')
        #pred_Y = regressor.predict(Test_X, num_iteration=regressor.best_iteration)
        pred_Y = regressor.predict(Test_X)
        # Save prediction and groud truth
        pred = np.array([Test_Y,pred_Y]).T
        
        cor = ma.corrcoef(ma.masked_invalid(pred_Y), ma.masked_invalid(Test_Y))[0,1]
        print("Prediction-gold standard's correlation on test set:",cor)
    else:
        print("Empty Test Set!")
        cor = np.nan
        pred = np.nan
    return cor, pred
        
def evaluate_SHAP_analysis(data_path, subset_type, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem):
    """ carry out subset-specific SHAP analysis
    
    Params
    ------
    data_path: str
        'test_by_cell_line' or 'test_by_indication'
    subset_type: str
        the subset division used to carry out SHAP:
        for example:
            '.'
            'moa': mode-of-action specific
            'tissue': tissue specific
            'cell_line': cell line specific
    target: str
        predicted score ('aoc' or 'bliss')

    dep_gene: str
        gene to make dependency plot
    
    Yields
    ------
    """
    from glob import glob

    all_path = sorted(list(glob(data_path+'/fold_*')))
    
    for path in all_path:
        
        idx = path.split('_')[-1]
        
        # load trained model parameter
        model_path = target+'_'+str(idx)+'.model'

        # define path of subsets
        if subset_type in ['Train', 'Test']:
            all_subset_path = glob(path+'/'+subset_type+'.tsv')
        else:
            all_subset_path = glob(path+'/'+subset_type+'/*.tsv')
        assert len(all_subset_path)>0, 'File path does not exist!'
    
        os.makedirs('./'+subset_type, exist_ok = True)

        for sub_path in all_subset_path:
            
            subset_name = sub_path.split('/')[-1].split('.')[0]

            if subset_type in ["Train", "Test"]:
                subset_outpath = './'+subset_type

            else:
                # all_moa           
                if subset_type == "all_moa":
                    subset_name = sub_path.split('/')[-1].split('.')[0]
                    if subset_name not in ['ATMi_ATRi', 'ATRi_PARPi', 'ATRi_TOP1i', 'ATRi_Cytostatic_Antimetabolite', 'DNAPKi_IR']:
                        continue
                    
                subset_outpath = './'+subset_type+'/'+subset_name
                os.makedirs(subset_outpath, exist_ok = True)
                
            print("Loading sub dataset for SHAP analysis: ",subset_name)
            Test = pd.read_csv(sub_path, sep = '\t')
            #print(Test.shape)
            #Test = Test[Test['.response_'+target].notna()]
            #print(Test.shape)
            if Test.shape[0] == 0:
                continue

            results = open(subset_outpath+'/cor_'+target+'.txt', 'a')
            # prediction on test set
            cor, pred = predict(sub_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path, surrogate_gene, surrogate_geneset, surrogate_chem)
            results.write('%.4f\n'%(cor))
            results.close()
            print("base value:", np.mean(pred[:,1]))

            # SHAP analysis on the test set
            shap_val, _, f_name = SHAP_analysis(sub_path,exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features,  model_path, surrogate_gene, surrogate_geneset, surrogate_chem)

            # save original Test_X without target gene and network altherations:
            if ('network' in features) or ('target_gene' in features):
                features = [i for i in features if (i not in ['network', 'target_gene'])]
            
            print(features)
            f_name, Test_X, Test_Y = build_feature_dataset(Test, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem, if_train = False)

            # Save SHAP values and feature values to np arrays
            shap_val = sparse.csr_matrix(shap_val)
            Test_X = sparse.csr_matrix(Test_X)
            sparse.save_npz(subset_outpath+'/SHAP_'+target+'_'+str(idx)+'.npz', shap_val)
            sparse.save_npz(subset_outpath+'/features_'+target+'_'+str(idx)+'.npz', Test_X)
            np.savetxt(subset_outpath+'/gs_'+target+'_'+str(idx)+'.txt', Test_Y)
            np.savetxt(subset_outpath+'/pred_'+target+'_'+str(idx)+'.txt', pred[:,1])
            pickle.dump(f_name,open(subset_outpath+'/feature_names_'+target+'.pkl','wb'))

def five_fold_cross_validation(data_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem):
    """ Five fold cross validation of synergy prediction model

    Params
    ------
    data_path: str
    exclude_synergy_batch: boolean
        use synergy batch as feature or not
    exclude_cell_line: boolean
    exclude_cancer_type: boolean
    target: str
        predicted score (aoc or bliss)
    features: a list of strings
        features used to construct feature set
    surrogate: for surrogate model
    
    Yields
    ------
    eva_df: Pandas dataframe
        evaluation results from k-fold cv
    eva_all: a str
        evalution reaulta from all k-fold; with 95 CI
    """
    from glob import glob
    
    eva_df = {"fold":[],"Pearson's r":[]} # predictions from all k-fold
    pred_all = [] # predictions from all k-fold
    
    all_path = sorted(list(glob(data_path+'/fold_*')))
    for path in all_path:

        idx = path.split('_')[-1]
        
        # start model training 
        train_path = path+'/Train.tsv'
        model_path = target+'_'+str(idx)+'.model'
        cor = train(train_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path, surrogate_gene, surrogate_geneset, surrogate_chem)

        # prediction on test set
        test_path = path+'/Test.tsv'
        cor, pred = predict(test_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path, surrogate_gene, surrogate_geneset, surrogate_chem)
        #np.savetxt(target+'_pred_'+str(idx)+'.txt', pred)
        
        # include fold prediction into all prediction
        pred_all.extend(pred.tolist())

        # write out test performance
        eva_df["fold"].append(idx)
        eva_df["Pearson's r"].append(cor)

    eva_df = pd.DataFrame.from_dict(eva_df)

    # calculate overall prediction performance confidence interval
    pred_all =  np.array(pred_all)
    mb, lb, ub = boostrapping_confidence_interval(pred_all, 0.95)
    eva_all = "mean[95CI]: %.4f[%.4f, %.4f]" % (mb, lb, ub)
    
    return eva_df, eva_all

def predict_on_held_out(data_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, surrogate_gene, surrogate_geneset, surrogate_chem):
    """ After model training, evaluate on the the held-out set.

    Params:
    """
    os.makedirs('hold_out', exist_ok = True)
    eva_df = {"fold":[],"Pearson's r":[]}
    pred_all = [] 
    for i in range(5):
        model_path = target+'_'+str(i)+'.model'
        # prediction on the held-out set
        cor, pred = predict(data_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path, surrogate_gene, surrogate_geneset, surrogate_chem)
        pred_all.extend(pred.tolist())
        
        print(cor)
        eva_df["fold"].append(i)
        eva_df["Pearson's r"].append(cor)
        np.savetxt('./hold_out/'+target+'_pred_'+str(i)+'.txt', pred)
    
    eva_df = pd.DataFrame.from_dict(eva_df)
    eva_df.to_csv('./hold_out/'+'cor_'+target+'.csv')
    pred_all =  np.array(pred_all)
    mb, lb, ub = boostrapping_confidence_interval(pred_all, 0.95)
    eva_all = "mean[95CI]: %.4f[%.4f, %.4f]" % (mb, lb, ub)

    result = open("./hold_out/cv-95ci_"+target+".txt", 'w')
    result.write(eva_all)
    result.close()

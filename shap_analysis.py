import pandas as pd
import numpy as np
from scipy.stats import pearsonr,spearmanr
import pickle
import json
from tqdm import tqdm
import lightgbm as lgb
import shap
import matplotlib.pyplot as plt
import glob
import os
from build_feature_dataset import *
from utils import *

def SHAP_analysis(test_path,exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, cv_split, features, mol_type, model_path, surrogate_gene, surrogate_geneset, surrogate_synleth, surrogate_chem):
    """ Carry out SHAP analysis of trained regressor on target dataset

    Parameters:
    -----------
    test_path:
    exclude_synergy_batch:
    target:
    features:
    model_path:

    Yields:
    -------
    shap_df: Pandas dataframe; n row by m column
        SHAP contribution of all features on all instances
    shap_fig: matplotlib plot file
    """
    
    print("Carrying out SHAP analysis on the test set...")
    regressor = pickle.load(open(model_path, 'rb'))

    Test = pd.read_csv(test_path, sep = '\t')
    #Test = Test[Test['.response_'+target].notna()]
    f_name, Test_X, Test_Y = build_feature_dataset(Test, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, cv_split, features, mol_type, surrogate_gene, surrogate_geneset, surrogate_synleth, surrogate_chem, if_train = False)
    
    # get feature values
    shap_values = shap.TreeExplainer(regressor).shap_values(Test_X)
    #shap_df = pd.DataFrame(np.array(shap_values), columns = f_name)
    
    # get direction of feature value on SHAP
    #cors = []
    #for i in range(shap_values.shape[1]):
    #    cor,_ = spearmanr(shap_values[:,i], Test_X[:,i])
    #    cors.append(cor)
    #shap_dir = {'features':f_name, 'cor':cors}

    # summary plot
    #shap.summary_plot(shap_values, Test_X, feature_names = f_name, show=False)
    #shap_fig = plt.gcf()
    #plt.close()

    #return shap_dir, shap_df, shap_fig
    return shap_values, Test_X, f_name

def Dependency_plot(gene, features, shap_values):
    """ Make dependency plot of different types of molecular markers (cnv, snv, exp, lof) of a single gene

    Parameteres:
    ------------
    gene: gene of interest
    
    Yields:
    -------
    dep_fig: dependency plot
    """
    main_marker = "exp"
    dependent_markers = ["cnv", "snv", "lof"]
    for dm in dependent_markers:
        main_f = gene+'_'+main_marker
        dep_f = gene+'_'+dm
        shap.dependence_plot(main_f, shap_values, features, interaction_index=dep_f)
        dep_fig = plt.gcf()
        plt.close()

    return dep_fig



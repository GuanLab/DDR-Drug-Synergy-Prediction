#!usr/bin/env python3
import pandas as pd
import numpy as np
from glob import glob
import pickle

from build_feature_dataset import *

def get_top_efficacy(data, model_path, target):
    """ calculate efficacy improved by pout predictions (compared to random assignmmnt)
    """
    exclude_synergy_batch = True
    exclude_cell_line = True
    exclude_cancer_type = True
    # for new cell lines
    # output predictions for all drug pairs available
    # get the best drug pair which ranks top
    # compare best efficacy to average efficacy

    # load model
    regressor = pickle.load(open(model_path, 'rb'))

    # load test set
    feature_set = ["monotherapy","moa","molecular","target_gene","geneset","chemical_structure","drug_name"]
    #feature_set = ["monotherapy","moa","molecular","target_gene","geneset","chemical_structure"]

    f_name, Test_X, Test_Y = build_feature_dataset(data, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features = feature_set, if_train = False)
    
    # get predictions
    pred_Y = regressor.predict(Test_X)
    # get top prediction index
    top_idx = np.argmax(pred_Y)
    # top_efficacy
    top_eff = Test_Y[top_idx]
    return top_eff

target = 'aoc'
all_path = glob.glob('../../test_by_cell_line/fold*')
improve_df = {'fold':[], 'cell_line':[], 'tissue':[], 'median_efficacy':[], 'selected_efficacy':[], 'improvement':[]}
all_improvement = []
for path in sorted(list(all_path)):
    idx = path.split('_')[-1]
    data_path = path+'/Test.tsv'
    model_path = 'monotherapy_'+target+'_'+str(idx)+'.model'
    data = pd.read_csv(data_path,sep ='\t')
    data = data[data['.response_'+target].notna()]
    all_cell_lines = data['.identifier_sample_name'].unique()
    for cell_line in all_cell_lines:
        sub_data = data.loc[data['.identifier_sample_name'] == cell_line,:]
        tissue = sub_data['.metadata_cancer_type'].iloc[0]
        print(tissue)
        ave_efficacy = sub_data['.response_'+target].median()
        top_efficacy = get_top_efficacy(sub_data, model_path, target)
        improvement = top_efficacy/ave_efficacy-1
        improvement_per = "{0:.0%}".format(improvement)
        improve_df['fold'].append(idx)
        improve_df['tissue'].append(tissue)
        improve_df['cell_line'].append(cell_line)
        improve_df['median_efficacy'].append(ave_efficacy)
        improve_df['selected_efficacy'].append(top_efficacy)
        improve_df['improvement'].append(improvement_per)
        if improvement>=0:
            all_improvement.append(1)
        else:
            all_improvement.append(0)

print("total_improvement",np.mean(np.array(all_improvement)))
improve_df = pd.DataFrame.from_dict(improve_df)
improve_df.to_csv("improvemt_efficacy.csv", index =False)

#!usr/bin/python3
#author @rayezh
import pandas as pd
from collections import defaultdict
feature_path = '../data/Monotherapy.tsv'
monotherapy = pd.read_csv(feature_path, sep = '\t' ,dtype = str)
"""
features:
    .identifier_batch
    .identifier_sample_name
    .metadata_moa_1
    .metadata_moa_2
    .metadata_treatment_remarks
    .metadata_treatment_1
    .metadata_treatment_2
    .response_aoc
    .response_bliss
    .metadata_cancer_type
    .metadata_cancer_subtype
"""
#clean unused columns
#assume all have 13 batches; for every batch, maximum 2 reps; moa as a feature 
a = defaultdict(lambda:np.zeros(27))
new_monotherapy = pd.DataFrame({'.identifier_batch':monotherapy['.identifier_batch'], 
    '.metadata_treatment_remarks': monotherapy['.metadata_treatment_remarks'],
    '.response_aoc':monotherapy['.response_aoc'],
    '.metadata_cancer_type':monotherapy['.metadata_cancer_type'],
    '.metadata_cancer_subtype':monotherapy['.metadata_cancer_subtype']})
moa = []
cell_drug = []
for _, row in monotherapy.iterrows():
    if row['.metadata_treatment_1'] == "Monotherapy":
        name = row['.metadata_treatment_2']
        moa.append(row['.metadata_moa_2'])
    else:
        name = row['.metadata_treatment_1']
        moa.append(row['.metadata_moa_1'])
    name = row['.identifier_sample_name']+'.'+name
    cell_drug.append(name)

new_monotherapy['.moa'] = moa
new_monotherapy['.cell_drug'] = cell_drug
        
new_monotherapy.to_csv('cleaned_mono.tsv', sep = '\t', index = False, na_rep = 'NA')

######################################## construct baseline model: baseline by cell_line-drug pair #########################################

new_monotherapy = pd.read_csv('cleaned_mono.tsv', sep = '\t')
cell_drug = list(set(list(new_monotherapy['.cell_drug'])))


"""
categorical features:
    .metadata_cancer_types
    .metadata_cancer_subtype
    .moa
    .metadata_treatment_remarks
pair: rep 26times
    .response_aoc(.identifier_batch)

"""

class embeddings:
    def __init__(self, data):

        names = list(set(data['.metadata_cancer_type']))
        self.cancer_type = {names[i]:i for i in range(len(names))}
        names = list(set(data['.metadata_cancer_subtype']))
        self.cancer_subtype =  {names[i]:i for i in range(len(names))}
        names = list(set(data['.moa']))
        self.moa = {names[i]:i for i in range(len(names))}
        names = list(set(data['.metadata_treatment_remarks']))
        self.remark = {names[i]:i for i in range(len(names))}

cat_embeddings = embeddings(new_monotherapy)

for i in cell_drug:
    features = open('./cell_drug_features/'+i, 'w')
    rows = new_monotherapy[new_monotherapy['.cell_drug'] == i] 
    rows.reset_index(drop = True, inplace = True)
    #print(cat_embeddings.cancer_type[rows['.metadata_cancer_type'][0]],cat_embeddings.cancer_subtype[rows['.metadata_cancer_subtype'][0]])
    features.write('%d\t' % cat_embeddings.cancer_type[rows['.metadata_cancer_type'][0]])
    features.write('%d\t' % cat_embeddings.cancer_subtype[rows['.metadata_cancer_subtype'][0]])
    features.write('%d\t' % cat_embeddings.moa[rows['.moa'][0]])
    features.write('%d\t' % cat_embeddings.remark[rows['.metadata_treatment_remarks'][0]])
    
    for j in range(13):
        batch_rows = rows[rows['.identifier_batch'] == 'Oncolead_'+'{:0>3}'.format(j)]
        batch_rows.reset_index(drop = True, inplace = True)
        r,c = batch_rows.shape
        i_rep = 2
        k = 0
        while k < r:
            features.write('%.4f\t' % batch_rows['.response_aoc'][k])
            k+=1
        while k < i_rep:
            features.write('%.4f\t' % 0)
            k+=1
    features.write('\n')
    features.close()





#construct nextstep model: feature by drug










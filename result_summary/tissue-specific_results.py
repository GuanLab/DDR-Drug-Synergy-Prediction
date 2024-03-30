import pandas as pd
import numpy as np
import numpy.ma as ma
from glob import glob
from utils import *

out_df = {'tissue':[], 'score':[], '95CI_mb':[], '95CI_lb':[], '95CI_ub':[]}
for score in ['aoc', 'bliss']:
    path = "../experiment_test_by_indication/features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure_drug_name/"+score+"_pred_*.txt"
    files = glob(path)
    for f in files:
        tissue = f.split('/')[-1].split('_')[2].split('.')[0]
        all_pred = np.loadtxt(f)
        pred_Y = all_pred[:,0]
        Test_Y = all_pred[:,1]
        #cor = ma.corrcoef(ma.masked_invalid(pred_Y), ma.masked_invalid(Test_Y))[0,1] # ignore nan
        mb, lb, ub = boostrapping_confidence_interval(all_pred, 0.95)
        #print(tissue, mb, lb, ub)
        out_df['tissue'].append(tissue)
        out_df['score'].append(score)
        out_df['95CI_mb'].append(mb)
        out_df['95CI_lb'].append(lb)
        out_df['95CI_ub'].append(ub)

out_df = pd.DataFrame(out_df)
# round to 4 decimal places
out_df = out_df.round(4)
out_df.to_csv('tissue_specific_results.csv', index=False)

#cervix,aoc,0.6611,0.6410,0.6818
#kidney,aoc,0.5454,0.5271,0.5739
#cervix,bliss,0.7822,0.7627,0.8064
#kidney,bliss,0.8429,0.8278,0.8585



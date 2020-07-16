import pandas as pd
from tqdm import tqdm
def sum_SHAP(target):
    df = {'feature':[], 'SHAP_val':[], 'rep':[], 'feature_type': []}
    for i in range(5):
        a = pd.read_csv('SHAP_'+target+'_'+str(i)+'.csv', header = 0)
        for col in tqdm(a.columns, total = len(a.columns.to_list())):
            f = col
            s = sum([abs(i) for i in a[col]])
            exp = str(i)
            df['feature'].append(f)
            df['SHAP_val'].append(s)
            df['rep'].append(exp)
            if col.startswith('Monotherapy1_.metadata_cancer') or col.startswith('Monotherapy2_.metadata_cancer'):
                df['feature_type'].append('Tissue type')
            elif col.startswith('Monotherapy1_batch') or col.startswith('Monotherapy2_batch'):
                df['feature_type'].append('Response')
            elif col.startswith('Monotherapy1_.moa') or col.startswith('Monotherapy2_.moa'):
                df['feature_type'].append('Mode-of-action')
            elif col.startswith('Monotherapy1_.metadata_treatment_remarks') or col.startswith('Monotherapy2_.metadata_treatment_remarks'):
                df['feature_type'].append('Remark')
            elif col.startswith('Monotherapy1_MACCS') or col.startswith('Monotherapy2_MACCS'):
                df['feature_type'].append('Chemical_Structure_MACCS')
            elif col.startswith('Monotherapy1_Morgan') or col.startswith('Monotherapy2_Morgan'):
                df['feature_type'].append('Chemical_Structure_Morgan')
            elif col.startswith('Monotherapy1_RDK') or col.startswith('Monotherapy2_RDK'):
                df['feature_type'].append('Chemical_Structure_RDK')
            elif col.startswith('Monotherapy1_FP2') or col.startswith('Monotherapy2_FP2'):
                df['feature_type'].append('Chemical_Structure_FP2')
            elif col.startswith('Monotherapy1_FP3') or col.startswith('Monotherapy2_FP3'):
                df['feature_type'].append('Chemical_Structure_FP3')
            elif col.startswith('Monotherapy1_FP4') or col.startswith('Monotherapy2_FP4'):
                df['feature_type'].append('Chemical_Structure_FP4')
            elif  col.endswith('_snv'):
                df['feature_type'].append('Molecular_snv')
            elif col.endswith('_cnv'):
                df['feature_type'].append('Molecular_cnv')
            elif col.endswith('_exp'):
                df['feature_type'].append('Molecular_exp')
            elif col.endswith('_lof'):
                df['feature_type'].append('Molecular_lof')
            elif col.endswith('__coh_pat'):
                df['feature_type'].append('Molecular_coh_pat')
            elif col.endswith('__lof_pat'):
                df['feature_type'].append('Molecular_lof_pat')
            elif col.endswith('_ddr'):
                df['feature_type'].append('Molecular_ddr')
            elif col == 'Synergy_batch':
                df['feature_type'].append('Synergy_batch')
            elif col.endswith('_drug'):
                df['feature_type'].append('Drug_category')
            else:
                df['feature_type'].append('Unknown')
    df = pd.DataFrame.from_dict(df)
    df.to_csv(target+'_SHAP.tsv',sep ='\t', index = False)

sum_SHAP(target = 'aoc')
sum_SHAP(target = 'bliss')

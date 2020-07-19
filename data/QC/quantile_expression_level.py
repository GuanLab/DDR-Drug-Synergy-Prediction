import pandas as pd
from glob import glob

def quantile_normalize_expression_levels(df):
    """
    Parameters:
    df: expression level feature 
    a pandas dataframe

    Yields:
    df_new: quantiled expression level feature
    a pandas dataframe
    """
    cols = df.columns
    exp_cols = [c for c in cols if c.endswith('_exp')]
    for idx, r in df.iterrows():
        m = r[exp_cols].mean()
        std = r[exp_cols].std()
        for c in exp_cols:
            df.loc[idx,c]= (r[c]-m)/std
    return df

path = '../merck_confidential_features_20200221.tsv'
data = pd.read_csv(path,sep = '\t',index_col=0 ,header = 0)
new_data = quantile_normalize_expression_levels(data)
new_data.to_csv('../quantiled_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
#fpaths = glob('../*_merck_confidential_features_20200221.tsv')
#for path in fpaths:
#    data = pd.read_csv(path,sep = '\t',index_col=0 ,header = 0)
#    if path.startswith('../exp'):
#        quantile_normalize_expression_levels(data)


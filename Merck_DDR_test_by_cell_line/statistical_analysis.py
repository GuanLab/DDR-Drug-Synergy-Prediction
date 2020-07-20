#!usr/bin/env/python3
import pandas as pd
import numpy as np
from collections import defaultdict
import json
from tqdm import tqdm
from glob import glob
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager

class SHAP_Analysis:
    def __init__(self, pred_target, moa_pair = None):
        """
        * 1. Summarize SHAP values of all features of machine learning models on specific mode-of-action dataset.
            This function will also generate a .tsv file for plotting figures in R.
        * 2. Subtract SHAP values of molecular markers of ddr genes.
        * 3. Highlight tissue-specific biomarkers for synergy prediction models
        
        -------------------------------------------------------------------------------------------
        Parameters:
            pred_target: prediction target ('aoc' or 'bliss')
                a string
            moa_pair: the mode-of-action(MOA) subset to analysis (e.g. 'ATMi_ATRi')
                if None, carry our analysis on overall dataset.
                None or a string

        Yields:
            self.pred_target: pred_target
            
            self.moa_pair: moa_pair

            self.fpaths: input path of SHAP analysis results from model_training.py
                a list of string
            self.cat_outpath: output path of categoricalized SHAP analysis results.
                a string
            self.mol_outpath: output path of molecular marker's SHAP analysis results.
                a string
            self.fig_outpath: output path of interaction heatmap of top features.
                a string
        """
        self.pred_target = pred_target
        self.moa_pair = moa_pair
        if self.moa_pair == None:
            self.fpaths = glob('SHAP_'+self.pred_target+'_*.csv')
            self.cat_outpath = self.pred_target+'_SHAP.tsv'
            self.mol_outpath = self.pred_target+'_molecular_SHAP.tsv'
            self.fig_outpath = self.pred_target+'_interaction_SHAP.pdf'
        else:
            self.fpaths = glob('./moa/'+moa_pair+'/'+'SHAP_'+pred_target+'_*.csv')
            self.cat_outpath = './moa/'+moa_pair+'/'+pred_target+'_SHAP.tsv'
            self.mol_outpath = './moa/'+moa_pair+'/'+pred_target+'_molecular_SHAP.tsv'
            self.fig_outpath = './moa/'+moa_pair+'/'+self.pred_target+'_interaction_SHAP.pdf'
        
    def Annote_feature_by_category(self):
        """ 
        * 1. Annotate SHAP feature set by feature categories.
            This function will generate a .tsv file for plotting figures in R.  
        * 2. Calculate interactions between feature and tissue ('cancer_type' and 'cancer_subtype')
                                                    
        Yields:
            self.category_df: Summarized SHAP contributions by category.
                a Pandas dataframe
                witten to <pred_target>+"_SHAP.tsv"
        """

        def feature_to_category(f):
            """
            * Annote feature by category.

            Parameters:
                f: name of input feature.
                    a string
    
            Yields:
                cat: category of the input feature.
                    a string
            """
            cat = 'Unknown'

            if f.startswith('Monotherapy1_.metadata_cancer') or f.startswith('Monotherapy2_.metadata_cancer'):
                cat = 'Tissue type'
            elif (f == 'Synergy_batch'):
                cat = 'Synergy_batch'
            elif f.startswith('Monotherapy1_batch') or f.startswith('Monotherapy2_batch'):
                cat = 'Monotherapy drug efficacy (AOC)'
            elif f.startswith('Monotherapy1_.moa') or f.startswith('Monotherapy2_.moa'):
                cat = 'Mode-of-action'
            elif (f == 'Monotherapy1_drug') or (f == 'Monotherapy2_drug'):
                cat = 'Drug name'
            elif f.startswith('Monotherapy1_.metadata_treatment_remarks') or f.startswith('Monotherapy2_.metadata_treatment_remarks'):
                cat = 'Concentration of Static Compound'
            elif f.startswith('Monotherapy1_MACCS') or f.startswith('Monotherapy2_MACCS'):
                cat = 'Chemical_Structure_MACCS'
            elif f.startswith('Monotherapy1_Morgan') or f.startswith('Monotherapy2_Morgan'):
                cat = 'Chemical_Structure_Morgan'
            elif f.startswith('Monotherapy1_RDK') or f.startswith('Monotherapy2_RDK'):
                cat = 'Chemical_Structure_RDK'
            elif f.startswith('Monotherapy1_FP2') or f.startswith('Monotherapy2_FP2'):
                cat = 'Chemical_Structure_FP2'
            elif f.startswith('Monotherapy1_FP3') or f.startswith('Monotherapy2_FP3'):
                cat = 'Chemical_Structure_FP3'
            elif f.startswith('Monotherapy1_FP4') or f.startswith('Monotherapy2_FP4'):
                cat = 'Chemical_Structure_FP4'
            elif  f.endswith('_snv'):
                cat = 'Molecular_snv'
            elif f.endswith('_cnv'):
                cat = 'Molecular_cnv'
            elif f.endswith('_exp'):
                cat = 'Molecular_exp'
            elif f.endswith('_lof'):
                cat = 'Molecular_lof'
            elif f.endswith('__coh_pat'):
                cat = 'Molecular_coh_pat'
            elif f.endswith('__lof_pat'):
                cat = 'Molecular_lof_pat'
            elif f.endswith('_ddr'):
                cat = 'Molecular_ddr'
            else:
                print(f)
    
            return cat

        
        print('Start prepare SHAP dataframe for '+self.pred_target+' score on overall dataset...')

        df = {'feature':[], 'SHAP_val':[], 'rep':[], 'feature_type': [], 'interation_with_tissue_type':[], 'interation_with_tissue_subtype':[]}
        for path in self.fpaths:
            data = pd.read_csv(path, header = 0)
            rep = str(path.split('_')[-1].split('.')[0]) #replcate index
            
            # calculate tissue correlation
            s_tissue_type = np.array([i for i in data['Monotherapy1_.metadata_cancer_type']])
            s_tissue_subtype = np.array([i for i in data['Monotherapy1_.metadata_cancer_subtype']])
            
            for col in tqdm(data.columns, total = len(data.columns.to_list())):
                feature = col #feature type
                shap_val = sum([abs(i) for i in data[col]]) #SHAP value: sum of absolue shape value across all samples
                shap_all = np.array([i for i in data[col]])
                cor_type, _ = pearsonr(shap_all, s_tissue_type)
                cor_subtype, _ =  pearsonr(shap_all, s_tissue_subtype)
                df['feature'].append(feature)
                df['SHAP_val'].append(shap_val)
                df['rep'].append(rep)
                cat = feature_to_category(feature)
                df['feature_type'].append(cat)
                df['interation_with_tissue_type'].append(cor_type)
                df['interation_with_tissue_subtype'].append(cor_subtype)
    
        self.category_df = pd.DataFrame.from_dict(df)
        self.category_df.to_csv(self.cat_outpath, sep = '\t', index = False)
        return self.category_df

    def Annote_molecular_features_of_genes(self, category_df):
        """ 
        * Subtract SHAP values of molecular markers of ddr genes.
    
        Yields:
            self.molecular_df: SHAP analysis of ddr genes' molecular markers
                a Pandas dataframe
                written to <pred_target>+"_molecular_SHAP.tsv"
        """
        print('Start tissue-specific analysis of all featurese for '+self.pred_target+' score prediction models ...')

        target_genes = open('../data/target_genes.txt', 'r').read().rstrip().split('\n')
        data_mol = category_df[category_df['feature_type'].isin(['Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof'])]
        df = {'gene': [],'mol_type':[], 'SHAP_val':[], 'rep':[], 'target':[]}

        for _,r in data_mol.iterrows():
            gene,mol_type = r['feature'].split('_')
            shap_val = r['SHAP_val']
            rep = r['rep']
            df['mol_type'].append(mol_type)
            df['SHAP_val'].append(shap_val)
            df['rep'].append(rep)
            if gene in target_genes:
                df['target'].append('True')
                df['gene'].append(gene+'*')
            else:
                df['target'].append('False')
                df['gene'].append(gene)
    
        self.molecular_df = pd.DataFrame.from_dict(df)
        self.molecular_df.to_csv(self.mol_outpath, sep = '\t', index = False)

    def analysis_top_feature_interation_network(self,n_top):
        """
        * Interaction map between top molecular features

        n_top: top n features to plot correlation heatmap
            a integer

        Yields:
            Correlation heatmaps of top n features

        """
        #combine SHAP Analysis results from all cross-validation models
        all_data = []
        for path in self.fpaths:
            data = pd.read_csv(path, header = 0)
            all_data.append(data) 
            #rep = str(path.split('_')[-1].split('.')[0]) #replcate index
        data = pd.concat(all_data)
        
        #reorder by the correlation with tissue type feature
        cor = {'feature':[], 'tissue_cor':[], 'abs_tissue_cor':[]}
        s_tissue_type = np.array([i for i in data['Monotherapy1_.metadata_cancer_type']])
        for col in data.columns:
            s = np.array([i for i in data[col]])
            cor_type, _ = pearsonr(s, s_tissue_type)
            cor['feature'].append(col)
            cor['tissue_cor'].append(cor_type)
            cor['abs_tissue_cor'].append(abs(cor_type))
        cor = pd.DataFrame.from_dict(cor)
        cor = cor.nlargest(n_top, 'abs_tissue_cor')
        cor = cor.sort_values(by=['tissue_cor'], ascending =False)
        sub_columns = cor['feature'].to_list() #top n features
            
        sub_data = data[sub_columns]
        cor_matrix = sub_data.corr()
        # output correlation matrix of top features 
        outpath = self.pred_target+'_'+str(n_top)+'_feature_interaction.tsv'
        cor_matrix.to_csv(outpath, sep = '\t', header = True, index = True)

        
        # plot correlation heatmap between features

        fig = plt.figure(dpi=300, figsize = (20,20))
        sns.set_style({'font.family':'serif', 'font.serif':'Helvetica'})
        p = sns.heatmap(cor_matrix, square = True,cmap = 'seismic', xticklabels=1, yticklabels=1,vmin = -1, vmax = 1, center = 0, cbar_kws={'label': "Pearson's corrlation", 'shrink':0.5})
        p.tick_params(labelsize = 8)
        #p.figure.axes[-1].yaxis.label.set_size(8)
        p.axes.set_title("Heatmap of Top "+str(n_top)+" Molecular Biomarkers", fontdict={'fontsize':20, 'fontweight': 'medium'})
        fig = p.get_figure()
        #fig.savefig(rep+'_'+self.fig_outpath)
        fig.savefig(self.fig_outpath)
        plt.clf()

        

def main():
        
    # 1. Yield dataframe of SHAP contribution by categories 
    # and
    # 2. Yield dataframe of SHAP contributions of ddr genes
    # and
    # 3. Check tissue-specific ddr gene features.
    df = SHAP_Analysis(pred_target = 'aoc').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'aoc').Annote_molecular_features_of_genes(df)
    df = SHAP_Analysis(pred_target = 'bliss').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'bliss').Annote_molecular_features_of_genes(df)

    df = SHAP_Analysis(pred_target = 'aoc', moa_pair = 'ATMi_ATRi').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'aoc', moa_pair = 'ATMi_ATRi').Annote_molecular_features_of_genes(df)
    df = SHAP_Analysis(pred_target = 'bliss', moa_pair = 'ATMi_ATRi').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'bliss', moa_pair = 'ATMi_ATRi').Annote_molecular_features_of_genes(df)

    df = SHAP_Analysis(pred_target = 'aoc', moa_pair = 'ATMi_PARPi').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'aoc', moa_pair = 'ATMi_PARPi').Annote_molecular_features_of_genes(df)
    df = SHAP_Analysis(pred_target = 'bliss', moa_pair = 'ATMi_PARPi').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'bliss', moa_pair = 'ATMi_PARPi').Annote_molecular_features_of_genes(df)

    df = SHAP_Analysis(pred_target = 'aoc', moa_pair = 'ATRi_PARPi').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'aoc', moa_pair = 'ATRi_PARPi').Annote_molecular_features_of_genes(df)
    df = SHAP_Analysis(pred_target = 'bliss', moa_pair = 'ATRi_PARPi').Annote_feature_by_category()
    SHAP_Analysis(pred_target = 'bliss', moa_pair = 'ATRi_PARPi').Annote_molecular_features_of_genes(df)

    # 4. Analysis interactions between top features.
    SHAP_Analysis(pred_target = 'aoc').analysis_top_feature_interation_network(n_top = 50)
    SHAP_Analysis(pred_target = 'bliss').analysis_top_feature_interation_network(n_top = 50)

if __name__ == '__main__':
    main()






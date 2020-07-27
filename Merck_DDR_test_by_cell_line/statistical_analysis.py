#!usr/bin/env/python3
import pandas as pd
import numpy as np
from collections import defaultdict
import json
import re
from tqdm import tqdm
from glob import glob
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager

class SHAP_Analysis:
    def __init__(self, pred_target, moa_pair = None, n_top = 50):
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

            n_top: number of top genes for interaction analysis
                an integer

        Yields:
            self.pred_target: pred_target
            
            self.moa_pair: moa_pair

            self.n_top:

            self.fpaths: input path of SHAP analysis results from model_training.py
                a list of string

            self.processed_fpaths: input path of SHAP analysis results, with SHAP value of every gene calculated.
                a list of string

            self.cat_outpath: output path of categoricalized SHAP analysis results.
                a string

            self.mol_outpath: output path of molecular marker's SHAP analysis results.
                a string

            self.top_gene_outpath: output path of top gene's SHAP analysis results.
                a string

            self.fig_outpath: output path of interaction heatmap of top features.
                a string

            self.category_df: categoricalized SHAP analysis results.
                a Pandas dataframe

            self.molecular_df: molecular marker's SHAP analysis results.
                a Pandas dataframe

            self.top_gene_df: top gene's SHAP analysis results.
                a Pandas dataframe
            self.top_genes: top genes
                a list of string
        """
        self.pred_target = pred_target
        self.moa_pair = moa_pair
        self.n_top = n_top
        if self.moa_pair == None:
            self.outpath = './'
        else:
            self.outpath = './moa/'+moa_pair+'/'
        
        self.SHAPpath = self.outpath+'SHAP_'+self.pred_target+'_*.csv'
            
        self.processed_SHAPpath = self.outpath+'processed_SHAP_'+self.pred_target+'_*.csv'
        
        assert len(glob(self.SHAPpath)) > 0, 'No SHAP Analysis results of machine learning models found!'

        if len(glob(self.processed_SHAPpath)) == 0: 
            print("Sum SHAP values of genes ...")
            self.Sum_SHAP_of_genes()
            
        self.fpaths = glob(self.SHAPpath)
        self.processed_fpaths = glob(self.processed_SHAPpath)
        self.cat_outpath = self.outpath+self.pred_target+'_SHAP.tsv'
        self.mol_outpath = self.outpath+self.pred_target+'_molecular_SHAP.tsv'
        self.top_gene_outpath = self.outpath+self.pred_target+'_top_gene.tsv'
        self.tissue_specificity_outpath = self.outpath+self.pred_target+'-SHAP-tissue-specificity_of_top_genes.tsv'
        self.fig_outpath = self.outpath+self.pred_target+'_interaction_SHAP.pdf'

        if len(glob(self.cat_outpath)) == 0:
            self.category_df = self.Annote_feature_by_category()
        else:
            self.category_df = pd.read_csv(self.cat_outpath, sep = '\t', header = 0)
        
        if len(glob(self.mol_outpath)) == 0:
            self.molecular_df = self.Annote_molecular_features_of_genes()
        else:
            self.molecular_df = pd.read_csv(self.mol_outpath, sep = '\t', header = 0)
        
        if len(glob(self.top_gene_outpath)) == 0:
            self.top_gene_df = self.Get_top_n_genes()
        else:
            self.top_gene_df = pd.read_csv(self.top_gene_outpath, sep = '\t', header = 0)

        self.top_genes = [re.sub('\*$', '', gene) for gene in self.top_gene_df['gene']]

        self.Analyse_tissue_specificity_of_top_features()
        self.Analyse_top_feature_interation_network()

    def feature_to_category(self,f):
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
        elif f.endswith('_all'):
            cat = 'Molecular_all'
        else:
            print("Unknown feature found: "+f)
    
        return cat

    def Sum_SHAP_of_genes(self):
        """
        * Preprocess SHAP value results from models of each fold.

        Parameters:
            self.fpaths

        Yields:
            processed SHAP values, saved to <processed>_fpath
        """
        def Sum_all_mols(df):
            """
            * Sum all molecular markers (cnv, snv, exp and lof) of this gene
            
            Parameters:
                df: SHAP value table
                    a Pandas dataframe
            Yields:
                df: processed SHAP value table with extra columns: <gene>_all_mol
                    a Pandas dataframe
            """
            genes = []
            all_col = df.columns
            for col in all_col:
                if self.feature_to_category(col) in ['Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof']:
                    genes.append(col.split('_')[0])
            
            genes = set(genes)
            for g in tqdm(genes):
                gene_col = [col for col in all_col if col.split('_')[0] == g]
                df[g+'_all'] = df[gene_col].sum(axis = 1)
            
            return df

        print('Combining molecular markers of genes ...')
        
        for path in self.fpaths:
            data = pd.read_csv(path, header = 0)
            data = Sum_all_mols(data)
            data.to_csv('processed_'+path, index = False)

    def Get_top_n_genes(self):
        """ 
        * Get n top genes
        
        Parameters:
            self.molecular_df
            self.n_top
            self.top_gene_outpath

        Yields:
            top_genes: top n genes selected and their mean SHAP values across all folds
                a Pandas dataframe
        """
        df = self.molecular_df.loc[self.molecular_df['mol_type'] == 'all'][['gene', 'SHAP_val']]
        top_gene_df = df.groupby('gene', as_index=False).agg('mean').nlargest(self.n_top, 'SHAP_val')
        top_gene_df.to_csv(self.top_gene_outpath, sep = '\t', index = False)
        return top_gene_df

    def Annote_feature_by_category(self):
        """ 
        * Annotate SHAP feature set by feature categories.
            This function will generate a .tsv file for plotting figures in R.  
        
        Parameters:
            self.fpaths

        Yields:
            category_df: Summarized SHAP contributions by category.
                a Pandas dataframe
                written to <pred_target>+"_SHAP.tsv"
        """

        print('Start prepare SHAP dataframe for '+self.pred_target+' score on overall dataset...')
        df = {'feature':[], 'SHAP_val':[], 'rep':[], 'feature_type': []} 
        for path in self.processed_fpaths:
            data = pd.read_csv(path, header = 0)
            rep = str(path.split('_')[-1].split('.')[0]) #replcate index
            
            for col in tqdm(data.columns, total = len(data.columns.to_list())):
                feature = col #feature type
                shap_val = sum([abs(i) for i in data[col]]) #SHAP value: sum of absolue shape value across all samples
                df['feature'].append(feature)
                df['SHAP_val'].append(shap_val)
                df['rep'].append(rep)
                cat = self.feature_to_category(feature)
                df['feature_type'].append(cat)

        category_df = pd.DataFrame.from_dict(df)
        category_df.to_csv(self.cat_outpath, sep = '\t', index = False)
        return category_df

    def Annote_molecular_features_of_genes(self):
        """ 
        * Subtract SHAP values of molecular markers of ddr genes.
        
        Parameters:
            self.pred_target
            self.category_df

        Yields:
            molecular_df: SHAP analysis of ddr genes' molecular markers
                a Pandas dataframe
                written to <pred_target>+"_molecular_SHAP.tsv"
        """
        print('Start tissue-specific analysis of all featurese for '+self.pred_target+' score prediction models ...')

        target_genes = open('../data/target_genes.txt', 'r').read().rstrip().split('\n')
        data_mol = self.category_df[self.category_df['feature_type'].isin(['Molecular_all', 'Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof'])]
        df = {'gene': [],'mol_type':[], 'SHAP_val':[], 'rep':[], 'target':[]}

        for _,r in data_mol.iterrows():
            gene,mol_type = r['feature'].split('_')
            shap_val = r['SHAP_val']
            rep = r['rep']
            df['mol_type'].append(mol_type)
            df['SHAP_val'].append(shap_val)
            df['rep'].append(rep)
            if gene in target_genes: # whether is a direct target of drug
                df['target'].append('True')
                df['gene'].append(gene+'*')
            else:
                df['target'].append('False')
                df['gene'].append(gene)
    
        molecular_df = pd.DataFrame.from_dict(df)
        molecular_df.to_csv(self.mol_outpath, sep = '\t', index = False)    
        return molecular_df

    def Analyse_tissue_specificity_of_top_features(self):
        """
        * Show tissue-specificity of top genes.
            1. Calculate interaction between top genes and tissue type features:
            (cancer_type and cancer_subtype)

            2. show which tissue type correlates well with which top genes
        
        Parameters:
            self.processed_fpaths
            self.top_genes
            self.pred_target

        Yields:
            df: tissue specificity based on SHAP value correlation
            
        """
        print('Start analysing tissue-specificity of top genes ...')
        df = {'gene':[], 'interation_with_tissue_type':[], 'interation_with_tissue_subtype':[]}
        all_data = []
        for path in self.processed_fpaths:
            data = pd.read_csv(path, header = 0)
            all_data.append(data)
        data = pd.concat(all_data)

        cols = [gene+'_all' for gene in self.top_genes]
        data_gene = data[cols]
        data_gene.rename(columns = {gene+'_all':gene for gene in self.top_genes}, inplace = True)

        # calculate correlation of SHAP values between tissue feature and top genes.
        s_tissue_type = np.array([i for i in data['Monotherapy1_.metadata_cancer_type']])
        s_tissue_subtype = np.array([i for i in data['Monotherapy1_.metadata_cancer_subtype']])

        for col in tqdm(data_gene.columns, total = len(data_gene.columns.to_list())):
            gene = col
            shap_all = np.array([i for i in data_gene[col]])
            cor_type, _ = pearsonr(shap_all, s_tissue_type)
            cor_subtype, _ =  pearsonr(shap_all, s_tissue_subtype)
            df['gene'].append(gene)
            df['interation_with_tissue_type'].append(cor_type)
            df['interation_with_tissue_subtype'].append(cor_subtype)
        
        df = pd.DataFrame.from_dict(df)
        df.to_csv(self.tissue_specificity_outpath, sep = '\t', index = False)

        #########################start from here: how to define correlation with each tissue type respectively?

    def Analyse_top_feature_interation_network(self):
        """
        * Draw interaction map between top genes

        Parameters:
            self.pred_target
            self.propcessed_fpaths
            self.category_df
            self.top_genes

        Yields:
            correlation heatmaps of top genes
        """

        #combine SHAP Analysis results from all cross-validation models
        all_data = []
        for path in self.processed_fpaths:
            data = pd.read_csv(path, header = 0)
            all_data.append(data) 
        
        data = pd.concat(all_data)
        cols = [gene+'_all' for gene in self.top_genes]
        data_gene = data[cols]
        data_gene.rename(columns = {gene+'_all':gene for gene in self.top_genes}, inplace = True)
    
        cor_matrix = data_gene.corr()
        # output correlation matrix of top features 
        outpath = self.pred_target+'_top_gene_interaction.tsv'
        cor_matrix.to_csv(outpath, sep = '\t', header = True, index = True)
        
        # plot correlation heatmap between features

        fig = plt.figure(dpi=300, figsize = (20,20))
        sns.set_style({'font.family':'serif', 'font.serif':'Helvetica'})
        p = sns.heatmap(cor_matrix, square = True,cmap = 'seismic', xticklabels=1, yticklabels=1,vmin = -1, vmax = 1, center = 0, cbar_kws={'label': "Pearson's corrlation", 'shrink':0.5})
        p.tick_params(labelsize = 8)
        #p.figure.axes[-1].yaxis.label.set_size(8)
        p.axes.set_title("Correlation Heatmap of Top Molecular Biomarkers", fontdict={'fontsize':20, 'fontweight': 'medium'})
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
    SHAP_Analysis(pred_target = 'aoc')
    #SHAP_Analysis(pred_target = 'bliss')

    """
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
    """

if __name__ == '__main__':
    main()






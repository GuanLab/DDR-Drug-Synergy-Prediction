#!usr/bin/env/python3
import pandas as pd
import numpy as np
from scipy import sparse
from collections import defaultdict
import json
import re
import pickle
import os
import shap
from tqdm import tqdm
from glob import glob
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager

class Downstream_SHAP_Analysis:
    def __init__(self, pred_target, moa_pair = None, subset = '', n_top = 20):
        """ Downstream analysis from SHAP
        * 1. Summarize SHAP values of all features of machine learning models on specific mode-of-action dataset.
            This function will also generate a .tsv file for plotting figures in R.
        * 2. Subtract SHAP values of molecular markers of ddr genes.
        * 3. Analysis tissue specificity of top genes
        * 4. Analysis the interactions between top genes 
        
        -------------------------------------------------------------------------------------------
        Params:
            pred_target: prediction target ('aoc' or 'bliss')
                a string

            moa_pair: the mode-of-action(MOA) subset to analysis (e.g. 'ATMi_ATRi')
                if None, carry our analysis on overall dataset.
                None or a string

            n_top: number of top genes for interaction analysis
                an integer
                default: 20

        Yields:
            self.pred_target: pred_target
            
            self.moa_pair: moa_pair

            self.n_top: top genes in *_top_gene.tsv.

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
        self.subset = subset
        self.outpath = '../'+subset+'/'+moa_pair+'/'

        # Check if SHAP analysis results are in the right place
        self.all_SHAPpath = glob(self.outpath+'SHAP_'+self.pred_target+'_*.npz')
        assert len(self.all_SHAPpath) > 0, 'No SHAP Analysis results of machine learning models found!'
        self.all_fpath =  glob(self.outpath+'features_'+self.pred_target+'_*.npz')
        assert len(self.all_fpath) > 0, 'No feature values of machine learning models found!'
        self.all_gspath =  glob(self.outpath+'gs_'+self.pred_target+'_*.txt')
        assert len(self.all_gspath) > 0, 'No feature values of machine learning models found!'

        self.fnames = pickle.load(open(self.outpath+'feature_names_'+self.pred_target+'.pkl', 'rb'))

        # combine all shap from cv
        print('Loading SHAP ...')
        all_shap = []
        for f in sorted(self.all_SHAPpath):
            n = sparse.load_npz(f)
            all_shap.append(n)
        self.all_shap = sparse.vstack(all_shap)
        # combine all features from cv
        print('Loading features ...')
        all_features = []
        for f in sorted(self.all_fpath):
            n = sparse.load_npz(f)
            all_features.append(n)
        self.all_features = sparse.vstack(all_features)
        all_gs =[]
        # combine all gs from cv
        print('Loading gold standards ...')
        for f in sorted(self.all_gspath): 
            n = np.loadtxt(f)
            all_gs.append(n)
        self.all_gs = np.concatenate(all_gs)

        # calculate direction of shap with features
        self.SHAP_summary_path = self.outpath+'SHAP_summary_'+self.pred_target+'.csv'
        self.SHAP_summary_plot_path =  self.outpath+'SHAP_summary_scatter_'+self.pred_target+'.pdf'
        print('Make summary plot ...')
        SHAP_all_df, shap_fig = self.calculate_SHAP_summary(self.all_shap, self.all_features,  self.all_gs, self.fnames)
        SHAP_all_df.to_csv(self.SHAP_summary_path, index = False)
        shap_fig.savefig(self.SHAP_summary_plot_path, format='pdf', dpi=1200, bbox_inches='tight')

        # dissect and sum shap for each category
        print('Dissect SHAP for each type ...')
        
        # summed by gene
        data_gene = self.Sum_all_mol_of_gene(self.all_shap, self.fnames)
        data_gene.to_csv(self.outpath+'gene_SHAP_'+self.pred_target+'.csv', index = False)
        # summed by category
        data_cat = self.Sum_category(self.all_shap, self.fnames)
        data_cat.to_csv(self.outpath+'cat_SHAP_'+self.pred_target+'.csv', index = False)
        # chemical-sum of each fingerprints
        data_chem = self.Sum_each_chemical(self.all_shap, self.fnames)
        data_chem.to_csv(self.outpath+'chem_SHAP_'+self.pred_target+'.csv', index = False)
        # molecular-sum of each type of molecular (4+2+1)
        data_mol = self.Sum_each_molecular(self.all_shap, self.fnames)
        data_mol.to_csv(self.outpath+'mol_SHAP_'+self.pred_target+'.csv', index = False)
        # synthetic lethality
        data_syn = self.Sum_each_sythetic_lethality(self.all_shap, self.fnames)
        data_syn.to_csv(self.outpath+'syn_SHAP_'+self.pred_target+'.csv', index = False)
        

        # locally dissect SHAP for each cancer type and cancer subtype
        sub_df_all = []
        data_gene_all = []
        data_cat_all = []
        data_chem_all = []
        data_mol_all = []
        data_syn_all = []
        for tissue, sub_shap, sub_features, sub_gs in self.generate_tissue_specific():
            #outpath = self.outpath+tissue+'/'
            #os.makedirs(outpath, exist_ok = True)
            #SHAP_summary_path = outpath+'SHAP_summary_'+self.pred_target+'.csv'
            #SHAP_summary_plot_path =  outpath+'SHAP_summary_scatter_'+self.pred_target+'.pdf'
            sub_df, sub_fig = self.calculate_SHAP_summary(sub_shap, sub_features, sub_gs, self.fnames)
            sub_df['tissue'] = tissue
            # summed by gene
            data_gene = self.Sum_all_mol_of_gene(sub_shap, self.fnames)
            data_gene['tissue'] = tissue
            #data_gene.to_csv(outpath+'gene_SHAP_'+self.pred_target+'.csv', index = False)
            # summed by category
            data_cat = self.Sum_category(sub_shap, self.fnames)
            data_cat['tissue'] = tissue
            #data_cat.to_csv(outpath+'cat_SHAP_'+self.pred_target+'.csv', index = False)
            # chemical-sum of each fingerprints
            data_chem = self.Sum_each_chemical(sub_shap, self.fnames)
            data_chem['tissue'] = tissue
            #data_chem.to_csv(outpath+'chem_SHAP_'+self.pred_target+'.csv', index = False)
            # molecular-sum of each type of molecular (4+2+1)
            data_mol = self.Sum_each_molecular(sub_shap, self.fnames)
            data_mol['tissue'] = tissue
            #data_mol.to_csv(outpath+'mol_SHAP_'+self.pred_target+'.csv', index = False)
            # synthetic lethality
            data_syn = self.Sum_each_sythetic_lethality(sub_shap, self.fnames)
            data_syn['tissue'] = tissue
            #data_syn.to_csv(outpath+'syn_SHAP_'+self.pred_target+'.csv', index = False)

            sub_df_all.append(sub_df)
            
            data_gene_all.append(data_gene)
            data_cat_all.append(data_cat)
            data_chem_all.append(data_chem)
            data_mol_all.append(data_mol)
            data_syn_all.append(data_syn)
            
        sub_df_all = pd.concat(sub_df_all)
        data_gene_all = pd.concat(data_gene_all)
        data_cat_all = pd.concat(data_cat_all)
        data_chem_all = pd.concat(data_chem_all)
        data_mol_all = pd.concat(data_mol_all)
        data_syn_all = pd.concat(data_syn_all)
        
        sub_df_all.to_csv(self.outpath+'tissue_summary_SHAP_'+self.pred_target+'.csv', index = False)
        data_gene_all.to_csv(self.outpath+'tissue_gene_SHAP_'+self.pred_target+'.csv', index = False)
        data_cat_all.to_csv(self.outpath+'tissue_cat_SHAP_'+self.pred_target+'.csv', index = False)
        data_chem_all.to_csv(self.outpath+'tissue_chem_SHAP_'+self.pred_target+'.csv', index = False)
        data_mol_all.to_csv(self.outpath+'tissue_mol_SHAP_'+self.pred_target+'.csv', index = False)
        data_syn_all.to_csv(self.outpath+'tissue_syn_SHAP_'+self.pred_target+'.csv', index = False)

        """
        self.Analyse_top_feature_interation_network()
        """

    def feature_to_category(self,f):
        """
        * Annote feature by category.

        Params:
            f: name of input feature.
                a string
    
        Yields:
            cat: category of the input feature.
                a string
        """

        cat = 'Unknown'

        if (f == 'Cell_line'):
            cat = 'Cell_line'
        elif (f == 'Synergy_batch'):
            cat = 'Synergy batch'
        elif (f == 'Cancer_type'):
            cat = 'Cancer_type'
        elif (f == 'Cancer_subtype'):
            cat = 'Cancer_subtype'
        elif (f == "Concentration"):
            cat = 'Concentration of Static Compound'
        elif re.match(r'Treatment_[1|2]_Oncolead_[0-9]{3}', f) is not None:
            cat = 'Monotherapy drug efficacy (AOC)'
        elif re.match(r'Treatment_[1|2]_ave', f) is not None:
            cat = 'Monotherapy drug efficacy (AOC)'
        elif re.match(r'Treatment_[1|2]_moa', f) is not None:
            cat = 'Mode-of-action'
        elif re.match(r'Treatment_[1|2]_name', f) is not None:
            cat = 'Drug name'
        elif re.match(r'Treatment_[1|2]_RDK_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_RDK'
        elif re.match(r'Treatment_[1|2]_MACCS_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_MACCS'
        elif re.match(r'Treatment_[1|2]_Morgan_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_Morgan'
        elif re.match(r'Treatment_[1|2]_FP2_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_FP2'
        elif re.match(r'Treatment_[1|2]_FP3_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_FP3'
        elif re.match(r'Treatment_[1|2]_FP4_[0-9]+', f) is not None:
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
        elif f.startswith('Geneset_'):
            cat = 'Geneset'
        elif f.startswith('CRISPR_cancer_dependency_'):
            cat = 'CRISPR Cancer Dependency'
        elif f.startswith('Synleth_'):
            cat = 'Synthetic lethality'
        else:
            print("Unknown feature found: "+f)
    
        return cat

    def calculate_SHAP_summary(self, shap_vals, features, gs, fnames):
        """
        * Calculate feature value's correlation with SHAP value
        Params:
        all_shap
        all_features
        all_gs
        fnames

        Yields:
        SHAP_dir: dict
        """
        SHAP_all_df = {'feature':fnames, 'mean|SHAP val|':[], 'dir(Pearsonr-SHAP)':[], 'dir(Pearsonr-gs)':[]}#, 'dir(p)':[]}
        s = shap_vals.toarray()
        f = features.toarray()
        for i,n in enumerate(tqdm(fnames)):
            #cor,p = spearmanr(s[:,i],f[:,i])
            #cor,p = pearsonr(s[:,i],f[:,i])
            x = s[:,i]
            y = f[:,i]
            nan = np.logical_or(np.isnan(x), np.isnan(y))
            cor1 = np.corrcoef(x[~nan],y[~nan])[0,1]
            x = gs
            nan = np.logical_or(np.isnan(x), np.isnan(y))
            cor2 = np.corrcoef(x[~nan],y[~nan])[0,1]   # groundtruth
            e = abs(s[:,i]).mean()
            SHAP_all_df['mean|SHAP val|'].append(e)
            #SHAP_all_df['dir(Spearmanr)'].append(cor)
            SHAP_all_df['dir(Pearsonr-SHAP)'].append(cor1)
            SHAP_all_df['dir(Pearsonr-gs)'].append(cor2)
            #SHAP_all_df['dir(p)'].append(p)
        
        SHAP_all_df = pd.DataFrame.from_dict(SHAP_all_df).sort_values(by = 'mean|SHAP val|', ascending=False)
        #SHAP_all_df.to_csv(summary_path, index = False)

        # summary plot
        shap.summary_plot(s, f, feature_names = fnames, show=False)
        shap_fig = plt.gcf()
        plt.close()
        #shap_fig.savefig(summary_plot_path, format='pdf', dpi=1200, bbox_inches='tight')

        return SHAP_all_df, shap_fig

    def Sum_all_mol_of_gene(self, shap_vals, fnames):
        """
        * Sum all molecular markers (cnv, snv, exp and lof) of this gene
            
        Params:
            shap: SHAP values
            fn: feature names
            Yields:
        df: processed SHAP value table with extra columns: <gene>
        a Pandas dataframe
        """
        target_genes = open('../../../feature/target_gene/target_genes.txt', 'r').read().rstrip().split('\n') #put '*' after target genes
        genes = []
        for n in fnames:
            if self.feature_to_category(n) in ['Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof']:
                genes.append(n.split('_')[0])
    
        df_new = {'gene':[], 'mean|SHAP val|':[]}
        genes = sorted(set(genes))
        for g in genes:
            if g in target_genes:
                g_new = g+"*"
            else:
                g_new = g
            gene_i = [i for i,f in enumerate(fnames) if f.split('_')[0] == g]
            df_new['gene'].append(g_new)
            df_new['mean|SHAP val|'].append(abs(shap_vals[:,gene_i].sum(axis =1)).mean())
            
        df_new = pd.DataFrame.from_dict(df_new).sort_values(by = 'mean|SHAP val|', ascending=False)
        return df_new
        
    def Sum_category(self,shap_vals,fnames):
        """ 
        * SHAP of each category
        Mode-of-action: 'Mode-of-action'
        Drug-name: 'Drug name'
        Monotherapy: 'Monotherapy drug efficacy (AoC)'
        Molecular Biomarkders: 'Molecular_*'
        Geneset Annotations: 'Geneset'
        Chemical Structure Fingerprints: 'Chemical_Structure_*'

        Params:
            shap: sparse matrix
            fnames: list of feature names

        Yields:
            df_new: SHAP by category
        """
        categories = {"Mode-of-action": r'Mode-of-action',
            "Drug name": r'Drug*',
            "Monotherapy": r'Monotherapy*',
            "Molecular Biomarkders": r'Molecular_*',
            "Geneset Annotations": r'Geneset',
            "Chemical Structure": r'Chemical_Structure_*',
            "Synthetic lethality": r'Synthetic lethality'}
        
        df_new = {'category':[], 'mean|SHAP val|':[]}
        for k,v in categories.items():
            cat_i = [i for i,col in enumerate(fnames) if re.match(v,self.feature_to_category(col)) is not None]
            df_new['category'].append(k)
            df_new['mean|SHAP val|'].append(abs(shap_vals[:,cat_i].sum(axis = 1)).mean())

        df_new = pd.DataFrame.from_dict(df_new).sort_values(by = 'mean|SHAP val|', ascending=False)
        return df_new

    def Sum_each_chemical(self, shap_vals, fnames):
        """ SHAP if each chemical strucure fingerprints
            6: MACCS, Morgan, RDK, FP2, FP3 FP4
        Params:

        Yields:
        """
        df_new = {'chem':[], 'mean|SHAP val|':[]}
        # SHAP of six types of fingerprints respectively
        for chem_fp in ['Chemical_Structure_RDK', 'Chemical_Structure_MACCS', 'Chemical_Structure_Morgan', 'Chemical_Structure_FP2', 'Chemical_Structure_FP3', 'Chemical_Structure_FP4']:
            chem_i = [i for i,col in enumerate(fnames) if self.feature_to_category(col) == chem_fp]
            df_new['chem'].append(chem_fp.split("_")[-1])
            df_new['mean|SHAP val|'].append(abs(shap_vals[:,chem_i].sum(axis = 1)).mean())
            
        df_new = pd.DataFrame.from_dict(df_new).sort_values(by = 'mean|SHAP val|', ascending=False)
        return df_new

    def Sum_each_molecular(self, shap_vals, fnames):
        """ 
        * SHAP of each mol biomarkers
            4+2+1
        """
        df_new = {'Molecular':[], 'mean|SHAP val|':[]}
        # Sum of each types of molecular biomarkers respectively
        for mol_marker in ['Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof', 'Molecular_coh_pat', 'Molecular_lof_pat', 'Molecular_ddr']:
            mol_i = [i for i,col in enumerate(fnames) if self.feature_to_category(col) == mol_marker]
            df_new['Molecular'].append("_".join(mol_marker.split("_")[1:]))
            df_new['mean|SHAP val|'].append(abs(shap_vals[:,mol_i].sum(axis = 1)).mean())

        df_new = pd.DataFrame.from_dict(df_new).sort_values(by = 'mean|SHAP val|', ascending=False)
        return df_new

    def Sum_each_sythetic_lethality(self, shap_vals, fnames):
        """
        """
        df_new = {"Synthetic_lethality":[], 'mean|SHAP val|':[]}
        syn_i = [i for i, col in enumerate(fnames) if self.feature_to_category(col) == "Synthetic lethality"]
        df_new["Synthetic_lethality"] = [fnames[i].split('Synleth_')[1] for i in syn_i]
        for i in syn_i:
            df_new['mean|SHAP val|'].append(abs(shap_vals[:,i]).mean())

        df_new = pd.DataFrame.from_dict(df_new).sort_values(by = 'mean|SHAP val|', ascending=False)
        return df_new

    def generate_tissue_specific(self):
        """
        * Locally analysis SHAP on each tissue (TODO: Cell lines)
        
        Params:
            self.gene_fpaths
            self.top_genes
            self.pred_target

        Yields:
            df: tissue specificity based on SHAP value correlation
            
        """
        # substract SHAP for each tissue (cancer type and cancer subtype)
        all_Test = []
        for f in sorted(self.all_SHAPpath):
            idx = f.split('_')[-1].split('.')[0]
            # look into data row idex 
            if self.moa_pair == '':
                testpath = '../../../test_by_cell_line/fold_'+idx+'/Test.tsv'   # this is for cell line model only! TODO: make this to cross-indication too
            else:
                testpath = '../../../test_by_cell_line/fold_'+idx+'/'+self.subset+'/'+self.moa_pair+'.tsv'
            print(testpath)
            Test = pd.read_csv(testpath, sep = '\t')
            all_Test.append(Test)
        all_Test = pd.concat(all_Test)
        all_Test = all_Test[all_Test['.response_'+self.pred_target].notna()].reset_index(drop = True)
        for tissue in sorted(set(all_Test['.metadata_cancer_type'])): # '.metadata_cancer_subtype'
            idx = all_Test.index[all_Test['.metadata_cancer_type'] == tissue]
            sub_shap = self.all_shap[idx,:]
            sub_features = self.all_features[idx,:]
            sub_gs = self.all_gs[idx]
            yield tissue, sub_shap, sub_features, sub_gs

    # TODO: synthetic lethality
    def Analyse_top_feature_interation_network(self):
        """
        * Draw interaction map between top genes

        Params:
            self.pred_target
            self.propcessed_fpaths
            self.category_df
            self.top_genes

        Yields:
            correlation heatmaps of top genes
        """

        #combine SHAP Analysis results from all cross-validation models
        all_data = []
        for path in glob(self.gene_SHAPpath):
            data = pd.read_csv(path, header = 0)
            all_data.append(data) 
        
        data = pd.concat(all_data)
        cols = list(self.top_genes)
        data_gene = data[cols]
        #data_gene.rename(columns = {gene+'_all':gene for gene in self.top_genes}, inplace = True)
    
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
        
    # global top features from SHAP 
    
    path = glob('../Test/')
    for p in path:
        tab = p.split('/')
        moa_pair = tab[2]
        subset = tab[1]
        Downstream_SHAP_Analysis(pred_target = 'aoc', moa_pair = moa_pair, subset= subset)
        Downstream_SHAP_Analysis(pred_target = 'bliss', moa_pair = moa_pair, subset= subset)

    path = glob('../Train/')
    for p in path:
        tab = p.split('/')
        moa_pair = tab[2]
        subset = tab[1]
        Downstream_SHAP_Analysis(pred_target = 'aoc', moa_pair = moa_pair, subset= subset)
        Downstream_SHAP_Analysis(pred_target = 'bliss', moa_pair = moa_pair, subset= subset)
    
    """
    #TO DO: fix these code as global
    # ATM/ATR/DNAPK-* combination treatment specific top features from SHAP
    moas = glob('../moa/*')
    for m in moas:
        moa = m.split('/')[-1]
        print(moa)
        Downstream_SHAP_Analysis(pred_target = 'aoc',  moa_pair = moa, subset= 'moa')
        Downstream_SHAP_Analysis(pred_target = 'bliss', moa_pair = moa, subset = 'moa')

    """
    
    # moa-pair treatment specific top features from SHAP
    moa_pairs = glob('../all_moa/*')
    for m in moa_pairs:
        tab = m.split('/')
        subset = tab[1]
        moa_pair = tab[2]
        if moa_pair in ['ATMi_ATRi', 'ATRi_PARPi', 'ATRi_TOP1i', 'ATRi_Cytostatic_Antimetabolite', 'DNAPKi_IR']:
            print(moa_pair)
            #Downstream_SHAP_Analysis(pred_target = 'aoc',  moa_pair = moa_pair, subset= subset)
            Downstream_SHAP_Analysis(pred_target = 'bliss',  moa_pair = moa_pair, subset= subset)


if __name__ == '__main__':
    main()






#!usr/bin/env/python3
#author: @rayezh
import pandas as pd
from glob import glob
import re
import numpy as np 

# TODO: five fold results
def cor_results_mols(split, score):
    # must include monotherapy
    cor_path = "../experiment_test_by_"+split+"/features_*only*/cor_"+score+".csv"
    allpath = glob(cor_path)
    all_data = []
    for path in allpath:
        mol_type = path.split("only")[-1].split("/")[0][1:]
        print(mol_type)
        data= pd.read_csv(path)
        data['features'] = mol_type
        all_data.append(data)
    all_data = pd.concat(all_data)
    all_data.to_csv("cor_all_mol"+split+"_"+score+".csv")

def cor_cv95_mols(split, score):
    cor_path = "../experiment_test_by_"+split+"/features_*only*/cv-95ci_"+score+".txt"
    allpath = glob(cor_path)
    df_all = {"feature":[], "mean[95CI]":[]}
    for path in allpath:
        mol_type = path.split("only")[-1].split("/")[0][1:]
        print(mol_type)
        p = re.sub(r'moa', 'mode-of-action', path) 
        f =  "+".join(p.split("/")[2].split("_")[1:])
        s = open(path, 'r').read().rstrip().split(":")[1]
        f = mol_type
        df_all["feature"].append(f)
        df_all["mean[95CI]"].append(s)
    df = pd.DataFrame.from_dict(df_all)
    df.to_csv("cor_meanci95_mol"+split+"_"+score+".csv")

def cor_results_chem(split, score):
    cor_path = "../experiment_test_by_"+split+"/features_chemical_structure_*/cor_"+score+".csv"
    allpath = glob(cor_path)
    all_data = []
    for path in allpath:
        mol_type = path.split("chemical_structure_")[-1].split("/")[0]
        print(mol_type)
        data= pd.read_csv(path)
        data['features'] = mol_type
        all_data.append(data)
    all_data = pd.concat(all_data)
    all_data.to_csv("cor_all_chem_"+split+"_"+score+".csv")

def cor_cv95_chem(split, score):
    cor_path = "../experiment_test_by_"+split+"/features_chemical_structure_*/cv-95ci_"+score+".txt"
    allpath = glob(cor_path)
    df_all = {"feature":[], "mean[95CI]":[]}
    for path in allpath:
        mol_type = path.split("chemical_structure_")[-1].split("/")[0]
        print(mol_type)
        s = open(path, 'r').read().rstrip().split(":")[1]
        f = mol_type
        df_all["feature"].append(f)
        df_all["mean[95CI]"].append(s)
    df = pd.DataFrame.from_dict(df_all)
    df.to_csv("cor_meanci95_chem_"+split+"_"+score+".csv")

def cor_results(split, score):
    # must include monotherapy
    #cor_path = "../experiment_test_by_"+split+"/features_geneset_surrogate_top*/cor_"+score+".csv"
    cor_path = "../experiment_test_by_"+split+"/features_*/cor_"+score+".csv"
    allpath = glob(cor_path)
    all_data = []
    for path in allpath:
        data= pd.read_csv(path)
        mol_type = path.split("top")[-1].split("/")[0][:]
        data['features'] = mol_type
        all_data.append(data)
    all_data = pd.concat(all_data)
    all_data.to_csv("cor_all_"+split+"_"+score+".csv")




def cor_cv95(split, score, if_holdout=False, if_topgenes=False, if_topsyn=False, if_only_mol=False):
    df_five_fold =  []
    df_all = {"features":[], "95CI_mb":[], '95CI_lb':[], '95CI_ub':[]}
    if if_holdout:
        test = "external validation"
        df_all = {"features":[], "95CI_mb":[], '95CI_lb':[], '95CI_ub':[], 'hold-out data':[]}
    else:
        test = "cross validation"
    
    f_path = "../experiment_test_by_"+split+"/features_*/"
    if if_topgenes:
        f_path = "../experiment_test_by_"+split+"/features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_topgene*/"
        df_all["topgene"] = []
        
    if if_topsyn:
        f_path = "../experiment_test_by_"+split+"/features_monotherapy_moa_synthetic_lethality_target_gene_network_geneset_chemical_structure_drug_name_topsynleth*/"
        df_all["topsyn"] = []
    
    if if_only_mol:
        f_path = "../experiment_test_by_"+split+"/features_molecular_target_gene_network_*/"
        df_all["mol_type"] = []

    allpath = sorted(glob(f_path))

    for path in allpath:
        f =  "+".join(path.split("/")[2].split("_")[1:])
        f = f.replace('target+gene', 'target_gene').replace('chemical+structure', 'chemical_structure').replace('synthetic+lethality', 'synthetic_lethality').replace('drug+name', 'drug_name').replace('coh+pat', 'coh_pat').replace('lof+pat', 'lof_pat')
        #f = path.split("top")[-1].split("/")[0][:]
        print(f)
        if if_topgenes:
            topgene = int(path.split("topgene")[1].split("/")[0])
        if if_topsyn:
            topsyn = int(path.split("topsynleth")[1].split("/")[0])
            print(topsyn)
        if if_only_mol:
            mol_type = f.split("+")[-1]
            print(mol_type)        
        # TODO: work on hold out data
        if if_holdout:
            holdout_path = path + "/hold_out/*"
            h_paths = glob(holdout_path)
            for h_path in h_paths:
                try:
                    cor_path = h_path+"/cv-95ci_"+score+".txt"
                    s = open(cor_path, 'r').read().rstrip().split(": ")[1]
                    cor = s.split('[')[0]
                    CI95 = s.split('[')[1].split(']')[0]
                    lb = CI95.split(',')[0]
                    ub = CI95.split(',')[1]
                    df_all["features"].append(f)
                    df_all["95CI_mb"].append(cor)
                    df_all["95CI_lb"].append(lb)
                    df_all["95CI_ub"].append(ub)
                    df_all["hold-out data"].append(h_path.split("/")[-1])
                    
                    # five fold
                    fold_path = h_path+"/cor_"+score+".csv"
                    data = pd.read_csv(fold_path)
                    data['features'] = f
                    data['hold-out data'] = h_path.split("/")[-1]
                    if if_topgenes:
                        df_all["topgene"].append(topgene)
                        data["topgene"] = topgene
                    if if_topsyn:
                        df_all["topsyn"].append(topsyn)
                        data["topsyn"] = topsyn
                    if if_only_mol:
                        df_all["mol_type"].append(mol_type)
                        data["mol_type"] = mol_type

                    df_five_fold.append(data)


                except:
                    pass             

        else:
            try:
                cor_path = path+"/cv-95ci_"+score+".txt"
                s = open(cor_path, 'r').read().rstrip().split(": ")[1]
                cor = s.split('[')[0]
                CI95 = s.split('[')[1].split(']')[0]
                lb = CI95.split(',')[0]
                ub = CI95.split(',')[1]
                df_all["features"].append(f)
                df_all["95CI_mb"].append(cor)
                df_all["95CI_lb"].append(lb)
                df_all["95CI_ub"].append(ub)

                # five fold
                fold_path = path+"/cor_"+score+".csv"
                data = pd.read_csv(fold_path)
                data['features'] = f
                if if_topgenes:    
                    df_all["topgene"].append(topgene)
                    data["topgene"] = topgene
                
                if if_topsyn:
                    df_all["topsyn"].append(topsyn)
                    data["topsyn"] = topsyn

                if if_only_mol:
                    df_all["mol_type"].append(mol_type)
                    data["mol_type"] = mol_type

                df_five_fold.append(data)
            except:
                pass
    
    print(df_all)
    df = pd.DataFrame.from_dict(df_all)
    # reorder df by correlation
    df = df.sort_values(by=['95CI_mb'], ascending=False).reset_index(drop=True)
    df['validation'] = test
    df['score'] = score
    df['split'] = split
    # set the highest performances to be *
    df['top'] = np.where(df['95CI_mb'] == df['95CI_mb'].max(), '*', '')
    #df.to_csv(outpath+"cor_meanci95_"+split+"_"+score+".csv")

    # five fold
    df_five_fold = pd.concat(df_five_fold)
    df_five_fold['validation'] = test
    df_five_fold['score'] = score
    df_five_fold['split'] = split

    if if_holdout:
        df_five_fold = df_five_fold.merge(df[['features', "95CI_mb", "95CI_lb", "95CI_ub", 'top', 'hold-out data']], on=['features', 'hold-out data'], how='left')
    else:
        df_five_fold = df_five_fold.merge(df[['features', "95CI_mb", "95CI_lb", "95CI_ub", 'top']], on=['features'], how='left')

    #df_five_fold.to_csv(outpath+"cor_five_fold_"+split+"_"+score+".csv")
    
    """
    from scipy.stats import ttest_rel
    for path1 in allpath:
        for path2 in allpath:
            fold_path1 = path1+"/cor_"+score+".csv"
            data1 = pd.read_csv(fold_path1)
            cor1 = data1["Pearson's r"]
            f1 = data1['features'][0]
            fold_path2 = path2+"/cor_"+score+".csv"
            data2 = pd.read_csv(fold_path2)
            cor2 = data2["Pearson's r"]
            f2 = data2['features'][0]
            # paired t-test
            t, p = ttest_rel(cor1, cor2, alternative='greater')
            if p < 0.05:
                print(f1, f2, p)
    """
    return df, df_five_fold


#df_all = []
#df_five_fold_all = []
#df_ho = []
#df_five_fold_ho= []
#df_topgene = []
#df_five_fold_topgene = []
#df_topgene_ho = []
#df_five_fold_topgene_ho = []
#df_five_fold_topsyn = []
#df_five_fold_topsyn_ho = []
df_five_fold_mol = []
df_five_fold_mol_ho = []

for i in ["cell_line","indication"]:  
    for j in ["aoc", "bliss"]:
        # all features
        #cor_results(i, j)
        #df1, df_five_fold1 = cor_cv95(i, j)
        #df_all.append(df1)
        #df_five_fold_all.append(df_five_fold1)
        # hold out
        #df2, df_five_fold2 = cor_cv95(i, j, if_holdout=True)
        #df_ho.append(df2)
        #df_five_fold_ho.append(df_five_fold2)
        # topgene only
        #df3, df_five_fold3 = cor_cv95(i, j, if_topgenes=True)
        #df_topgene.append(df3)
        #df_five_fold_topgene.append(df_five_fold3)
        # topgene only hold out
        #df4, df_five_fold4 = cor_cv95(i, j, if_topgenes=True, if_holdout=True)
        #df_topgene_ho.append(df4)
        #df_five_fold_topgene_ho.append(df_five_fold4)
        # topsyn only
        #df5, df_five_fold5 = cor_cv95(i, j, if_topsyn=True)
        #df_five_fold_topsyn.append(df_five_fold5)
        # topsyn only hold out
        #df6, df_five_fold6 = cor_cv95(i, j, if_topsyn=True, if_holdout=True)
        #df_five_fold_topsyn_ho.append(df_five_fold6)
        # molecular feature only
        df7, df_five_fold7 = cor_cv95(i, j, if_only_mol=True)
        df_five_fold_mol.append(df_five_fold7)
        # molecular feature only hold out
        df8, df_five_fold8 = cor_cv95(i, j, if_only_mol=True, if_holdout=True)
        df_five_fold_mol_ho.append(df_five_fold8)

#df_all = pd.concat(df_all)
#df_five_fold_all = pd.concat(df_five_fold_all)
#df_all.to_csv("cor_meanci95_all.csv")
#df_five_fold_all.to_csv("cor_five_fold_all.csv")

#df_ho = pd.concat(df_ho)
#df_five_fold_ho = pd.concat(df_five_fold_ho)
#df_ho.to_csv("cor_meanci95_ho.csv")
#df_five_fold_ho.to_csv("cor_five_fold_ho.csv")

#df_topgene = pd.concat(df_topgene)
#df_five_fold_topgene = pd.concat(df_five_fold_topgene)
#df_topgene.to_csv("cor_meanci95_topgene.csv")
#df_five_fold_topgene.to_csv("cor_five_fold_topgene.csv")

#df_topgene_ho = pd.concat(df_topgene_ho)
#df_five_fold_topgene_ho = pd.concat(df_five_fold_topgene_ho)
#df_topgene_ho.to_csv("cor_meanci95_topgene_ho.csv")
#df_five_fold_topgene_ho.to_csv("cor_five_fold_topgene_ho.csv")

#df_five_fold_topsyn = pd.concat(df_five_fold_topsyn)
#df_five_fold_topsyn.to_csv("cor_five_fold_topsyn.csv")
#df_five_fold_topgene_ho = pd.concat(df_five_fold_topsyn_ho)
#df_five_fold_topgene_ho.to_csv("cor_five_fold_topsyn_ho.csv")
        
df_five_fold_mol = pd.concat(df_five_fold_mol)
df_five_fold_mol.to_csv("cor_five_fold_mol.csv")
df_five_fold_mol_ho = pd.concat(df_five_fold_mol_ho)
df_five_fold_mol_ho.to_csv("cor_five_fold_mol_ho.csv")

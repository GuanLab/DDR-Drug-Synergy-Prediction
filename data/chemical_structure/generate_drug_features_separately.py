#!usr/bin/python3
import pandas as pd
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdmolops
from pubchempy import *
import openbabel
import pybel

def fp_to_feature_v2(fp):
    feature = [0 for i in range(1024)]
    for n in fp.bits:
        feature[n] = 1
    return feature
def fp_to_feature_v3(fp, max_bit):
    feature = [0 for i in range(max_bit+1)]
    for n in fp.bits:
        feature[n] = 1
    return feature

df_map = pd.read_csv('./merck_confidential_molecular_structures_20200420/smiles/smiles_and_ms_no.csv', header = 0) #treatment_name_orig --> treatment
max_bit_fp3 = 55
max_bit_fp4 = 307

for _,r in df_map.iterrows():
    smile = r['treatment_smiles']
    ms = Chem.MolFromSmiles(smile)
    mol = pybel.readstring("smi", smile)
    #MACCS features (167*1)
    fp = MACCSkeys.GenMACCSKeys(ms)
    tmp = fp.ToBitString()
    feature_1 = list(map(int, tmp))
    np.savetxt('./MACCS_features/'+r['treatment_name_orig'],np.array(feature_1))

    #Morgan Fingerprints(1024*1)
    fp = AllChem.GetMorganFingerprintAsBitVect(ms,2,nBits=1024)
    tmp = fp.ToBitString()
    feature_2 = list(map(int, tmp))
    np.savetxt('./Morgan_features/'+r['treatment_name_orig'],np.array(feature_2))

    #FP2 (1024*1)
    fp = mol.calcfp('FP2')
    feature_3 = fp_to_feature_v2(fp)
    np.savetxt('./FP2_features/'+r['treatment_name_orig'],np.array(feature_3))
    
    #FP3 (55*1)
    fp = mol.calcfp('FP3')
    feature_4 = fp_to_feature_v3(fp, max_bit_fp3)
    np.savetxt('./FP3_features/'+r['treatment_name_orig'],np.array(feature_4))

    #FP4 (307*1)
    fp = mol.calcfp('FP4')
    feature_5 = fp_to_feature_v3(fp, max_bit_fp4)
    np.savetxt('./FP4_features/'+r['treatment_name_orig'],np.array(feature_5))

    #RDK Fingerprints (2048*1)
    fp = rdmolops.RDKFingerprint(ms)
    tmp = fp.ToBitString()
    feature_6 = list(map(int, tmp))
    np.savetxt('./RDK_features/'+r['treatment_name_orig'],np.array(feature_6))

    the_feature = np.array(feature_1+feature_2+feature_3+feature_4+feature_5+feature_6)
    #print(the_feature.shape) #4627
    np.savetxt('./all_drug_features/'+r['treatment_name_orig'],the_feature)


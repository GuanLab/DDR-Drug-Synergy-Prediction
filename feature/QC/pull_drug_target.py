#!usr/env/python3
# pull drug targets from databases
import pandas as pd
import re
# ChEMBL
from chembl_webresource_client.new_client import new_client

#uniprot
import urllib.parse
import urllib.request

import json
import itertools

# generate ./target_genes

# all drug-response
df = pd.read_csv('../target_gene/all_drugs_summary.csv')
# find all drugs
all_drugs = df['drug_name']

# find all drugs with mode of action
drug2moa = {r['drug_name']:r['mode-of-action'] for _,r in df.iterrows()}

# all ddr genes
all_ddr_genes = open('../target_gene/ddr_gene.txt').read().strip().split('\n')
#print(all_ddr_genes)

comp2target = dict()

# LINCS drug interaction dataset:
# wget https://lincs.hms.harvard.edu/db/datasets/20000/results?search=&output_type=.csv
df_LINCS = pd.read_csv('./LINCS/dataset_20000_20210415165646.csv')

try:
    print("Start new compound-target list ...")
    comp2target = dict()
    for d in all_drugs:
        comp2target[d] = list()
except:
    pass
def uniprot2genesymbol(query):
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
            'from': 'ACC+ID',
            'to': '	GENENAME',
            'format': 'tab',
            'query': query
            }
    
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    response = response.decode('utf-8').strip().split('\n')[1:]
    symbol =[i.split('\t')[1] for i in response]
    return symbol

def target_from_moa(drug):
    """
    """
    genes = []
    v = drug2moa[drug]
    for v_sub in v.split('_'):
        genes.extend(re.findall("[a-zA-Z0-9]+i$", v_sub))
    genes = [i[:-1] for i in genes]
    moa_symbol = []
    for i in genes:
        r = re.compile(i+"[0-9]*$")
        new_genes = list(filter(r.match, all_ddr_genes))  # match genes in all_ddr_gene list
        # for the following genes, the protein (target name) and the gene names are different
        if i == 'DNAPK':
            new_genes = ['PRKDC']
        if i == 'CHK1':
            new_genes = ['CHEK1']
        if i == 'CHK2':
            new_genes = ['CHEK2']
        moa_symbol.extend(new_genes)

    return moa_symbol

def target_from_chembl(drug):
    """
    """
    molecule = new_client.molecule
    all_mols = []
    all_targets = []
    uniprots = set()
    for i in drug.split('_'):
        res = molecule.search(i)
        all_mols.extend(res)
    all_chembl_id = [i['molecule_chembl_id'] for i in all_mols]
    activities = new_client.target.filter(molecule_chembl_id__in=all_chembl_id).only(['target_chembl_id'])
    target_ids = [i['target_chembl_id'] for i in activities]
    targets = new_client.target.filter(target_chembl_id__in=target_ids).only(['target_components'])
    uniprots |= {comp['accession'] for t in targets for comp in t['target_components']}
    query = ' '.join(uniprots)
    symbol = uniprot2genesymbol(query)
    ddr_symbol = [i for i in symbol if i in all_ddr_genes]
    return ddr_symbol

def target_from_lincs(drug):
    """
    """
    # LINCS drug interaction dataset:
    # wget https://lincs.hms.harvard.edu/db/datasets/20000/results?search=&output_type=.csv
    lincs_symbol = df_LINCS[(df_LINCS['HMSL Small Mol Name'].str.contains(drug.replace('_','|'), case = False))&(df_LINCS['Binding Class'].isin([1,2,3]))]['HUGO Gene Symbol'].to_list()
    ddr_lincs_symbol = [i for i in lincs_symbol if i in all_ddr_genes]
    return ddr_lincs_symbol

def target_from_dgi(drug):
    """
    """
    # DGIdb: https://dgidb.org/api/v2/interactions.json?drugs=VELIPARIB
    dgi_symbol = []
    url = 'https://dgidb.org/api/v2/interactions.json?drugs='+drug.replace('_',',')
    with urllib.request.urlopen(url) as response:
        data = json.loads(response.read().decode())

    for i in range(len(data['matchedTerms'])):
            data_i = data['matchedTerms'][i]["interactions"] # interaction
            for j in range(len(data_i)):
                data_j = data_i[j]
                #print(data_j['interactionTypes'])
                if 'inhibitor' in data_j['interactionTypes']:
                    dgi_symbol.append(data_j['geneName'])
                #print(data_j['interactionTypes'], data_j['geneName'], d, drug2moa[d])
    ddr_dgi_symbol = [i for i in dgi_symbol if i in all_ddr_genes]
    return ddr_dgi_symbol

refs = []
for d in all_drugs:
    print(d, drug2moa[d])
    # find targets from MoA
    moa_symbol = target_from_moa(d)
    comp2target[d].extend(moa_symbol)
    ref = []
    """
    # find targets from ChEMBL
    ddr_symbol = target_from_chembl(d)
    comp2target[d].extend(ddr_symbol)
    if len(ddr_symbol) > 0:
        ref.append('CHEMBL')
    """
    # find targets from LINCS
    ddr_lincs_symbol = target_from_lincs(d)
    print('LINCS', ddr_lincs_symbol)
    comp2target[d].extend(ddr_lincs_symbol)
    if len(ddr_lincs_symbol) > 0:
        ref.append('LINCS')

    # find targets from DGIdb
    ddr_dgi_symbol = target_from_dgi(d)
    print('DGIdb', ddr_dgi_symbol)
    comp2target[d].extend(ddr_dgi_symbol)
    if len(ddr_dgi_symbol) > 0:
         ref.append('DGIdb')
    
    comp2target[d] = sorted(set(comp2target[d]))
    print(d,drug2moa[d],comp2target[d])
    ref = ';'.join(ref)
    refs.append(ref)


out_df = {'drug_name':list(comp2target.keys()),
        'mode-of-action':[drug2moa[i] for i in comp2target.keys()],
        'target_proteins':[','.join(i) for i in comp2target.values()],
        'references': refs}
out_df = pd.DataFrame.from_dict(out_df)
out_df.to_csv('../target_gene/all_drugs_summary.csv', index = False)

json.dump(comp2target,open('../target_gene/drug2gene.json', 'w'))
print(comp2target)
target_genes = '\n'.join(set(list(itertools.chain(*comp2target.values()))))
f = open('../target_gene/target_genes.txt', 'w')
f.write(target_genes)
f.close()

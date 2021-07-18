import os
import pandas as pd
import traversalHelper as tr
import numpy as np
from collections import defaultdict
from statistics import mean
import sys
pd.options.mode.chained_assignment = None

# data standard: nodeA, nodeB, type, score

def parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.176.tab2.txt', root='./'
    , wFile_GGI='./data/parsed/BioGRID_GGI.pkl', wFile_PPI='./data/parsed/BioGRID_PPI.pkl'):
    wFile_GGI, wFile_PPI, filename = root+wFile_GGI, root+wFile_PPI, root+filename

    if os.path.exists(wFile_GGI) and os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_GGI), pd.read_pickle(wFile_PPI)

    bioGRID_df = pd.read_csv(filename, sep='\t')
    colMap = {"Official Symbol Interactor A": "nodeA",
     "Official Symbol Interactor B": 'nodeB',
      "Experimental System Type": "type", "Score": "score"}
    
    ggi_df = bioGRID_df[bioGRID_df["Experimental System Type"] == 'genetic']
    ppi_df = bioGRID_df[bioGRID_df["Experimental System Type"] == 'physical']

    dfs = {'ggi': ggi_df, 'ppi': ppi_df}

    for key in dfs:
        for col in dfs[key].columns:
            if col not in colMap: dfs[key].drop(col, axis=1, inplace=True)
        dfs[key].rename(colMap, axis=1, inplace=True)
        dfs[key].reset_index(inplace=True)

    for key in dfs:
        sortedBR = tr.Helper.list_to_pathStrs([sorted(br) for br in np.asarray(dfs[key][['nodeA', 'nodeB']])])
        scoreMap, typeMap = defaultdict(list), {}
        for i in range(len(sortedBR)): # remove directional interaction
            if dfs[key]['nodeA'][i] == dfs[key]['nodeB'][i]: continue # remove self-interaction
            scoreMap[sortedBR[i]].append(np.nan if dfs[key]['score'][i] == '-' else float(dfs[key]['score'][i]))
            typeMap[sortedBR[i]] = dfs[key]['type'][i]
        
        for k, arr in scoreMap.items():
            if np.isnan(np.asarray(scoreMap[k])).all(): scoreMap[k] = np.nan
            else: scoreMap[k] = np.nanmean(np.asarray(scoreMap[k])) 

        ppi = list(np.transpose(np.asarray(tr.Helper.pathStrs_to_list(list(scoreMap.keys())))))
        dfs[key] = pd.DataFrame({'nodeA': ppi[0], 'nodeB': ppi[1],
         'type': list(typeMap.values()), 'score': list(scoreMap.values())})

    pd.to_pickle(dfs['ggi'], wFile_GGI)
    pd.to_pickle(dfs['ppi'], wFile_PPI)
    return dfs['ggi'], dfs['ppi']

def parse_geneName_map(filename="./data/BioGRID/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.176.tab2.txt", root="./"):
    #  systematic name => gene name
    bioGRID_df = pd.read_csv(root+filename, sep="\t")
    sysName = list(bioGRID_df['Systematic Name Interactor A'])+list(bioGRID_df['Systematic Name Interactor B'])
    geneName = list(bioGRID_df['Official Symbol Interactor A'])+list(bioGRID_df['Official Symbol Interactor B'])
    geneMap = dict(zip(geneName, sysName))
    reverseGeneMap = dict(zip(sysName, geneName))
    return geneMap, reverseGeneMap

if __name__ == "__main__":
    ggi_df, ppi_df = parse_bioGRID()
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction
    ppi = [sorted(arr) for arr in np.asarray(ggi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction

    ggi_df, ppi_df = parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.187.tab2.txt', root='./'
    , wFile_GGI='./data/parsed/BioGRID_homo_GGI.pkl', wFile_PPI='./data/parsed/BioGRID_homo_PPI.pkl')
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction
    ppi = [sorted(arr) for arr in np.asarray(ggi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction

import pandas as pd
import numpy as np
import os
import traversalHelper as tr
import sys

def parse_uniProt_map(uniProtMap):
    df = pd.read_csv(uniProtMap, sep='\t')
    df.dropna(inplace=True)
    uniProtMapping = dict(zip(list(df['Entry']), list(df['Gene names  (primary )'])))
    return uniProtMapping

def parse_MINT(ppiFile='./data/MINT/species yeast', uniProtMap="./data/UniProt/uniprot-taxonomy_559292.tab", 
                wFile_GGI='./data/parsed/MINT_GGI.pkl', wFile_PPI='./data/parsed/MINT_PPI.pkl', root='./'):
    ppiFile, uniProtMap, wFile_GGI, wFile_PPI = root+ppiFile, root+uniProtMap, root+wFile_GGI, root+wFile_PPI
    if os.path.exists(wFile_GGI) and os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_GGI), pd.read_pickle(wFile_PPI)

    uniProtMapping = parse_uniProt_map(uniProtMap)
    # only direct interaction & physical association & association are PPIs
    df = pd.read_csv(ppiFile, sep='\t', header=None)

    ppi_df = df[df[11].isin([i for i in df[11] if 'direct interaction' in i or 'physical association' in i or 'association' in i])].copy()
    ggi_df = df[~df[11].isin([i for i in df[11] if 'direct interaction' in i or 'physical association' in i or 'association' in i])].copy()
    dfs = {'ppi': ppi_df, 'ggi': ggi_df}
    
    for key in dfs:
        inv_ppis = np.transpose(np.asarray([dfs[key][0], dfs[key][1]]))
        ppis = []
        for i in inv_ppis:
            if i[0] == "-" or i[1] == "-": continue # remove empty entry
            if i[0] == i[1]: continue # remove self-interaction
            if not 'uniprotkb' in i[0] or not 'uniprotkb' in i[1]: continue #remove invalid entry (cannot map uniprot)
            ppi_h = i[0].split(":")[1].split("-")[0]
            ppi_t = i[1].split(":")[1].split("-")[0]
            ppis.append(sorted([ppi_h, ppi_t]))
        ppis = tr.Helper.pathStrs_to_list(set(
            tr.Helper.list_to_pathStrs(ppis))) # remove directional

        mappedPPIs = []
        for ppi in ppis:
            if ppi[0] not in uniProtMapping or ppi[1] not in uniProtMapping: continue
            if uniProtMapping[ppi[0]] == uniProtMapping[ppi[1]]: continue
            mappedPPIs.append(sorted([uniProtMapping[ppi[0]], uniProtMapping[ppi[1]]]))
        mappedPPIs = tr.Helper.pathStrs_to_list(set(
            tr.Helper.list_to_pathStrs(mappedPPIs)))
        mappedPPIs = np.transpose(np.asarray(mappedPPIs))
        dfs[key] = pd.DataFrame({'nodeA': mappedPPIs[0], 'nodeB': mappedPPIs[1]})

    pd.to_pickle(dfs['ppi'], wFile_PPI)
    pd.to_pickle(dfs['ggi'], wFile_GGI)
    return dfs['ggi'], dfs['ppi']


if __name__ == "__main__":
    ggi_df, ppi_df = parse_MINT()
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction
    ppi = [sorted(arr) for arr in np.asarray(ggi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction


    ggi_df, ppi_df = parse_MINT(ppiFile='./data/MINT/species human', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab", 
                wFile_GGI='./data/parsed/MINT_homo_GGI.pkl', wFile_PPI='./data/parsed/MINT_homo_PPI.pkl', root='./')
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction
    ppi = [sorted(arr) for arr in np.asarray(ggi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction

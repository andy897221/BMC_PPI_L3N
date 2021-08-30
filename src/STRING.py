import pandas as pd
import bioGRID as bg
import traversalHelper as tr
import numpy as np
import os
from collections import defaultdict
from statistics import mean
from scipy import stats

def parse_uniProt_map(uniProtMapF):
    df = pd.read_csv(uniProtMapF, sep='\t')
    df.dropna(inplace=True)
    uniProtMapping = dict(zip([i.split(";")[0] for i in df['Cross-reference (STRING)']], list(df['Gene names  (primary )'])))
    return uniProtMapping

def parse_STRING(ppiFile='./data/STRING/4932.protein.links.v11.0.txt'
    , typeFile='./data/STRING/4932.protein.actions.v11.0.txt'
    , uniProtMap='./data/UniProt/uniprot-taxonomy_559292_STRING.tab', root='./'
    , wFile_GGI='./data/parsed/STRING_GGI.pkl', wFile_PPI='./data/parsed/STRING_PPI.pkl'):
    ppiFile, typeFile, wFile_GGI, wFile_PPI, uniProtMap = root+ppiFile, root+typeFile, root+wFile_GGI, root+wFile_PPI, root+uniProtMap

    if os.path.exists(wFile_GGI) and os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_GGI), pd.read_pickle(wFile_PPI)

    # Sys name (used by STRING) => gene name (used by this project)
    reverseGeneMap = parse_uniProt_map(uniProtMap)

    df_STRING = pd.read_csv(ppiFile, sep=' ')
    dfType_STRING = pd.read_csv(typeFile, sep='\t')
    dfType_STRING.rename(columns={'item_id_a': 'protein1', 'item_id_b': 'protein2'}, inplace=True)
    merged_df = pd.merge(df_STRING, dfType_STRING, on=['protein1', 'protein2'])
    merged_df = merged_df[['protein1', 'protein2', 'combined_score', 'mode']]
    # merged_df['combined_score'] = merged_df['combined_score']/1000
    merged_df['combined_score'] = merged_df['combined_score']
    
    ppis = [list(arr) for arr in np.asarray(merged_df[['protein1', 'protein2']])]
    scoreList, typeList = list(merged_df['combined_score']), list(merged_df['mode'])
    rebuild_score, rebuild_type = defaultdict(list), defaultdict(list)
    for i in range(len(ppis)):
        ppi = ppis[i]
        if ppi[0] not in reverseGeneMap or ppi[1] not in reverseGeneMap: continue # convert to gene name
        ppi = sorted([reverseGeneMap[ppi[0]], reverseGeneMap[ppi[1]]]) # sort for remove directional
        if ppi[0] == ppi[1]: continue # remove self-interaction
        rebuild_score["\t".join(ppi)].append(scoreList[i])
        rebuild_type["\t".join(ppi)].append(typeList[i])
    
    # build ppi & ggi
    ppi_df, ggi_df = defaultdict(list), defaultdict(list)
    for ppi in rebuild_score:
        rebuild_score[ppi], rebuild_type[ppi] = np.asarray(rebuild_score[ppi]), np.asarray(rebuild_type[ppi])
        ppi_score, ggi_score = rebuild_score[ppi][rebuild_type[ppi] == 'binding'], rebuild_score[ppi][rebuild_type[ppi] != 'binding']
        ppi_type, ggi_type = rebuild_type[ppi][rebuild_type[ppi] == 'binding'], rebuild_type[ppi][rebuild_type[ppi] != 'binding']
        ppi = ppi.split("\t")
        if len(ppi_score) > 0:
            ppi_df['nodeA'].append(ppi[0])
            ppi_df['nodeB'].append(ppi[1])
            ppi_df['type'].append('binding')
            ppi_df['score'].append(np.mean(ppi_score))
        if len(ggi_score) > 0:
            ggi_df['nodeA'].append(ppi[0])
            ggi_df['nodeB'].append(ppi[1])
            ggi_df['type'].append(stats.mode(ggi_type)[0][0])
            ggi_df['score'].append(np.mean(ggi_score))

    ppi_df, ggi_df = pd.DataFrame(ppi_df), pd.DataFrame(ggi_df)
    pd.to_pickle(ggi_df, wFile_GGI)
    pd.to_pickle(ppi_df, wFile_PPI)
    return ggi_df, ppi_df

if __name__ == "__main__":
    ggi_df, ppi_df = parse_STRING()
    # print(set(list(STRING_GGI['type'])))
    # {'catalysis', 'ptmod', 'inhibition', 'activation', 'reaction', 'expression'}
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction
    ppi = [sorted(arr) for arr in np.asarray(ggi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction


    ggi_df, ppi_df = parse_STRING(ppiFile='./data/STRING/9606.protein.links.v11.0.txt'
    , typeFile='./data/STRING/9606.protein.actions.v11.0.txt'
    , uniProtMap='./data/UniProt/uniprot-organism__Homo+sapiens+[9606]_.tab', root='./'
    , wFile_GGI='./data/parsed/STRING_homo_GGI.pkl', wFile_PPI='./data/parsed/STRING_homo_PPI.pkl')
    print(ppi_df.head())
    print(ggi_df.head())
    print(len(ggi_df.index), len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction
    ppi = [sorted(arr) for arr in np.asarray(ggi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction

import pandas as pd
import numpy as np
import os
import traversalHelper as tr

def parse_uniProt_map(uniProtMap):
    df = pd.read_csv(uniProtMap, sep='\t')
    df.dropna(inplace=True)
    uniProtMapping = dict(zip(list(df['Entry']), list(df['Gene names  (primary )'])))
    return uniProtMapping

def parse_HuRI(ppiFile='./data/atlas/HuRI.psi', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab", 
                wFile_PPI='./data/parsed/HuRI_PPI.pkl', root='./'):
    ppiFile, uniProtMap, wFile_PPI = root+ppiFile, root+uniProtMap, root+wFile_PPI
    if os.path.exists(wFile_PPI): return pd.read_csv(wFile_PPI)

    uniProtMapping = parse_uniProt_map(uniProtMap)
    # only direct interaction & physical association & association are PPIs
    ppi_df = pd.read_csv(ppiFile, sep='\t', header=None)
    inv_ppis = np.transpose(np.asarray([ppi_df[0], ppi_df[1]]))
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
    ppi_df = pd.DataFrame({'nodeA': mappedPPIs[0], 'nodeB': mappedPPIs[1]})

    ppi_df.to_csv(wFile_PPI, index=False)
    return ppi_df


if __name__ == "__main__":
    ppi_df = parse_HuRI()
    print(ppi_df.head())
    print(len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction

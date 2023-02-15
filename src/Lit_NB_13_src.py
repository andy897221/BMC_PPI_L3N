import pandas as pd
import numpy as np
import os
import traversalHelper as tr

def parse_Lit_NB_src(ppiFile='./data/other github repo/L3-master/lit-nb-13.txt', 
wFile_PPI='./data/parsed/Lit_NB_src.pkl', root='./'):
    ppiFile, wFile_PPI = root+ppiFile, root+wFile_PPI
    if os.path.exists(wFile_PPI):
        return pd.read_pickle(wFile_PPI)

    ppi_df = pd.read_csv(ppiFile, sep="\t", header=None)
    inv_ppis = np.transpose(np.asarray([ppi_df[0], ppi_df[1]]))
    inv_ppis = [sorted(arr) for arr in inv_ppis]
    ppis = tr.Helper.pathStrs_to_list(set(
        tr.Helper.list_to_pathStrs(inv_ppis)
    ))
    ppis = np.transpose(np.asarray(ppis))
    ppi_df = pd.DataFrame({'nodeA': ppis[0], 'nodeB': ppis[1]})
    pd.to_pickle(ppi_df, wFile_PPI)
    return ppi_df

if __name__ == "__main__":
    ppi_df = parse_Lit_NB_src()
    print(ppi_df.head())
    print(len(ppi_df.index))
    ppi = [sorted(arr) for arr in np.asarray(ppi_df[['nodeA', 'nodeB']])]
    print(len(set(tr.Helper.list_to_pathStrs(ppi)))
    , len(tr.Helper.list_to_pathStrs(ppi)), len([p for p in ppi if p[0] == p[1]])) # check duplicate & self-interaction

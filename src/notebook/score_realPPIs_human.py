import sys, os
sys.path.append('..')

# standard
import json
from collections import defaultdict
import pandas as pd
import numpy as np
import seaborn as sns
from statistics import median
from matplotlib import pyplot as plt
from sklearn import metrics
import matplotlib.ticker as ticker
import random as rn

# my lib
import PPILinkPred as pred
import genData_helper as helper
import traversalHelper as tr
import bioGRID
import STRING
import MINT
import HuRI

ds_names = ['bioGRID_human', 'STRING_human', 'MINT_human', 'HuRI']
import_funcs = [
    bioGRID.parse_bioGRID(filename='./data/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.187.tab2.txt'
        , wFile_GGI='./data/parsed/BioGRID_homo_GGI.pkl'
        , wFile_PPI='./data/parsed/BioGRID_homo_PPI.pkl', root="../")

    , STRING.parse_STRING(ppiFile='./data/STRING/9606.protein.links.v11.0.txt'
        , typeFile='./data/STRING/9606.protein.actions.v11.0.txt'
        , uniProtMap='./data/UniProt/uniprot-taxonomy_9606_STRING.tab', root='../'
        , wFile_GGI='./data/parsed/STRING_homo_GGI.pkl', wFile_PPI='./data/parsed/STRING_homo_PPI.pkl')

    , MINT.parse_MINT(ppiFile='./data/MINT/species human', uniProtMap="./data/UniProt/uniprot-taxonomy_9606.tab"
        , wFile_GGI='./data/parsed/MINT_homo_GGI.pkl', wFile_PPI='./data/parsed/MINT_homo_PPI.pkl', root="../")
]
completePPIs_map = [
    [list(ppi) for ppi in np.asarray([*import_funcs[0]][1][['nodeA', 'nodeB']])]
    , [list(ppi) for ppi in np.asarray([*import_funcs[1]][1][['nodeA', 'nodeB']])]
    , [list(ppi) for ppi in np.asarray([*import_funcs[2]][1][['nodeA', 'nodeB']])]
    , [list(ppi) for ppi in np.asarray(HuRI.parse_HuRI(root="../")[['nodeA', 'nodeB']])]
]
ppi_ds = dict(zip(ds_names, completePPIs_map))

if __name__ == "__main__":
    methods = ['commonNeighbor', 'L3Normalizing', "CH2_L3", "Sim", "CRA"]+["L3E1_f{}".format(i) for i in range(1,3)]
    for ds in ppi_ds:
        for method in methods:
            predPPI, predScore = pred.multiCore_PPILinkPred(ppi_ds[ds], method, coreNo=14, logging=True
                                                            , testmode=True, nodePairs=ppi_ds[ds])
            print(ds, method, 'finished')
            with open('{}/{}_score.json'.format("ppiScoring_human_out", method+"_"+ds), 'w') as f:
                f.write(json.dumps(predScore))
            with open('{}/{}_ppi.json'.format("ppiScoring_human_out", method+"_"+ds), 'w') as f:
                f.write(json.dumps(predPPI))
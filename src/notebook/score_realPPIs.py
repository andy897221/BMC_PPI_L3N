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

# my lib
import PPILinkPred as pred
import genData_helper as helper
import traversalHelper as tr
import bioGRID
import STRING
import MINT

ppi_ds = {}
import_funcs = [bioGRID.parse_bioGRID(root="../"), STRING.parse_STRING(root="../"), MINT.parse_MINT(root="../")]
names = ['bioGRID', 'STRING', 'MINT']

for n in range(len(names)):
    _, df = import_funcs[n]
    ppi = [list(arr) for arr in np.asarray(df[['nodeA', 'nodeB']])]
    ppi_ds[names[n]] = ppi

if __name__ == "__main__":
    methods = ["commonNeighbor", "L3Normalizing", "CRA", "CH2_L3", "Sim"]+[j+"_"+i for j in ['L3E1'] for i in ['f1', 'f2']]
    for ds in ppi_ds:
        for method in methods:
            if os.path.exists("./ppiScoring_out/{}_{}_ppi.json".format(method, ds)): continue

            predPPI, predScore = pred.multiCore_PPILinkPred(ppi_ds[ds], method, coreNo=14, logging=True
                                                            , testmode=True, nodePairs=ppi_ds[ds])
            print(ds, method, 'finished')
            with open('{}/{}_score.json'.format("ppiScoring_out", method+"_"+ds), 'w') as f:
                f.write(json.dumps(predScore))
            with open('{}/{}_ppi.json'.format("ppiScoring_out", method+"_"+ds), 'w') as f:
                f.write(json.dumps(predPPI))
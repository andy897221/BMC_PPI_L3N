import sys
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
    methods = ['commonNeighbor', 'L3Normalizing', "CH2_L3", "Sim", "CRA"]
    methods += ["L3E{}_f{}".format(j, i) for j in range(1,5) for i in range(1,3)]
    for ds in ppi_ds:
        sampled_nonPPIs = [] # sampled in "generate prediction.ipynb"
        with open("sampled_datasets/{}_sampled_nonPPIs.json".format(ds), "r") as f:
            sampled_nonPPIs = json.loads(f.read())
        
        for trial in range(len(sampled_nonPPIs)):
            for method in methods:
                predPPI, predScore = pred.multiCore_PPILinkPred(
                    ppi_ds[ds], method, coreNo=14, logging=True, nodePairs=sampled_nonPPIs[trial])
                    
                print(ds, method, 'finished')
                with open('{}/{}_nonPPI_score_{}.json'.format("ppiScoring_out", method+"_"+ds, trial), 'w') as f:
                    f.write(json.dumps(predScore))
                with open('{}/{}_nonPPI_ppi_{}.json'.format("ppiScoring_out", method+"_"+ds, trial), 'w') as f:
                    f.write(json.dumps(predPPI))
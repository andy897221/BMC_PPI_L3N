# std
import numpy as np
import random as rn
import json
import time
from collections import defaultdict
from itertools import combinations

# datasets
import STRING
import MINT
import IntAct
import bioGRID
import HuRI

# my lib
import PPILinkPred as pred
import genData_helper as helper

tags = ['xyContrib_dualCN_uvJoin']
methods = ['interStr']
scoreArgDicts = [{'xyContrib': 'basic', 'dualCN': 'basic', 'uvJoin': 'basic'}]

# read dataset
samplePPIs = []
with open("./sampled_datasets/STRING_sampledPPIs.json", "r") as f:
    samplePPIs = json.loads(f.read())

for t in range(len(tags)):
    for i in range(len(samplePPIs)):
        saveFilename = 'ExactL31_STRING_sample_{}'.format(i)
        startTime = time.time()
        
        scoreArgsDict = scoreArgDicts[t]
        normOrder = ['normFunc', 'uvSpec', 'xySpec', 'uvContrib', 'xyContrib', 'dualCN', 'uvJoin']
        scoreArgs = ['null' if normTag not in scoreArgsDict else scoreArgsDict[normTag] for normTag in normOrder]

        samplePPIr = pred.ns.BRToRelat(pred.ns.toDualBR(samplePPIs[i]), rSet=True)
        sampleNodes = pred.ns.BRToNode(samplePPIs[i])
        nodePairs = list(combinations(sampleNodes, 2))

        pred._PPILinkPred(nodePairs, samplePPIr, "interStr", scoreArgs=scoreArgs, logging=False)
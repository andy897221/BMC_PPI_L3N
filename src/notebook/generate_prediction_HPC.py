import sys
sys.path.append('..')

# std
import numpy as np
import random as rn
import json
import time
import math
from collections import defaultdict
from itertools import combinations
from mpi4py import MPI

# datasets
import STRING
import MINT
import bioGRID
import HuRI

# my lib
import PPILinkPred as pred
import genData_helper as helper

# using HPC cores
def multiCore_PPILinkPred_HPC(samplePPIbr, scoringMethod, coreNo, topNo=None, logging=False, nodePairs=None):

    samplePPIr = pred.ns.BRToRelat(pred.ns.toDualBR(samplePPIbr), rSet=True)
    sampleNodes = pred.ns.BRToNode(samplePPIbr)
    if nodePairs is None: nodePairs = list(combinations(sampleNodes, 2))

    splitStartIndex = [i*math.floor(len(nodePairs)/coreNo) for i in range(0, coreNo)] # both splitting is correct
    splitEndIndex = [(i+1)*math.floor(len(nodePairs)/coreNo) if i != coreNo-1 else len(nodePairs) for i in range(0, coreNo)]
    if logging: logging = [True if i == 0 else False for i in range(coreNo)]
    else: logging = [False for i in range(coreNo)]

    rank = MPI.COMM_WORLD.Get_rank()
    nodePairs = nodePairs[splitStartIndex[rank]:splitEndIndex[rank]]
    scores, predictedPPIbrs = pred._PPILinkPred(nodePairs, samplePPIr, scoringMethod, logging)
    # if logging: print("\n")

    return predictedPPIbrs, scores, rank

def main():
    methods = ["L3E1_f1"] # examples, all data is generated part-by-part
    ds_names = ['bioGRID_human']

    for ds_name in ds_names:
        # read dataset
        samplePPIs = []
        with open("./sampled_datasets/{}_sampledPPIs.json".format(ds_name), "r") as f:
            samplePPIs = json.loads(f.read())

        # do link prediction & save results
        for method in methods:
            # for i in range(len(samplePPIs)):
            for i in range(0, 10):
                saveFilename = "{}_{}_sample_{}".format(method, ds_name, i)
                startTime = time.time()
                
                # jupyter notebook cannot display multi core logging, do it only in terminal
                predPPI, predScore, rank = multiCore_PPILinkPred_HPC(samplePPIs[i], method, coreNo=24, logging=False)
                if rank == 0: helper.write_runTime(saveFilename, time.time()-startTime)
                helper.write_resultData(predPPI, predScore, saveFilename+"_c{}".format(rank))

if __name__ == "__main__":
    main()
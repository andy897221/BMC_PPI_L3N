import traversalHelper as tr
from itertools import combinations
import helper as hr
import math, time
import multiprocessing
from multiprocessing import Pool, Manager
from functools import partial
import numpy as np
from statistics import mean
import json, sys, os
from collections import defaultdict
import random as rn

class helperFunc:
    @staticmethod
    def logging(count, lastCount, total, avgTime, startTime, frequency=1000):
        count += 1
        if count == 1: avgTime = time.time()-startTime
        else: avgTime = (avgTime*(count-1)+(time.time()-startTime))/count
        if count-lastCount > frequency:
            print("reference core's count: {}/{}. tick rate: {}. Expected sec to finish (hr): {} ({})".format(
                count, total, frequency, round(avgTime*(total-count), 2), round(avgTime*(total-count)/60/60, 2)), end="\r")
            lastCount = count
        return count, lastCount, avgTime


class ns:
    BRToRelat = tr.Helper.binary_to_relation
    toDualBR = tr.Helper.to_dual_binary_relation
    BRToNode = tr.Helper.binary_relation_to_node
    arr_pStr = tr.Helper.list_to_pathStrs
    pStr_arr = tr.Helper.pathStrs_to_list
    br_str = tr.Helper.br_to_pathStr
    L3Based = ["L3Normalizing", "L3Raw", "Sim"]
    L2Based = ["commonNeighbor"]
    CARBased = ["CRA", "CAR", "CH2_L3"]
    controlBased = ["countP4"]

class control:
    @staticmethod
    def handler(samplePPIr, nodeX, nodeY, scoringMethod):
        score = 0
        if scoringMethod == "countP4":
            uvPair, _, _ = L3.get_uv(nodeX, nodeY, samplePPIr)
            score = len(uvPair)
        
        return score

class L3:
    @staticmethod
    def L3_normalization(PPIr, uvPair, normFunc):
        score = 0
        for uv in uvPair:
            if normFunc(len(PPIr[uv[0]])*len(PPIr[uv[1]])) == 0: continue
            score += 1/normFunc(len(PPIr[uv[0]])*len(PPIr[uv[1]]))
        return score

    @staticmethod
    def Sim(samplePPIr, nodeX, nodeY, uvPair):
        nodeUs, nodeVs = set([uv[0] for uv in uvPair]), set([uv[1] for uv in uvPair])
        score = 0
        N = samplePPIr
        for v in nodeVs:
            score += (len(N[v]&N[nodeX])/len(N[v]|N[nodeX]))
        for u in nodeUs:
            score += (len(N[u]&N[nodeY])/len(N[u]|N[nodeY]))
        return score

    @staticmethod
    def get_uv(x, y, PPIr, uvJoin=True):
        candidateUs = PPIr[x]
        candidateVs = PPIr[y]
        if not uvJoin:
            candidateUs = candidateUs-candidateVs
            candidateVs = candidateVs-candidateUs
        uvPair = []
        for u in candidateUs:
            for v in candidateVs:
                if u not in PPIr[v]: continue
                uvPair.append([u,v])
        return uvPair, candidateUs, candidateVs

    @staticmethod
    def handler(samplePPIr, nodeX, nodeY, scoringMethod):
        if scoringMethod == "L3Normalizing":
            uvPair, candidateUs, candidateVs = L3.get_uv(nodeX, nodeY, samplePPIr)
            score = L3.L3_normalization(samplePPIr, uvPair, math.sqrt)
        elif scoringMethod == "L3Raw":
            uvPair, candidateUs, candidateVs = L3.get_uv(nodeX, nodeY, samplePPIr)
            score = L3.L3_normalization(samplePPIr, uvPair, lambda x: 1)
        elif scoringMethod == "Sim":
            uvPair, candidateUs, candidateVs = L3.get_uv(nodeX, nodeY, samplePPIr)
            score = L3.Sim(samplePPIr, nodeX, nodeY, uvPair)
        return score

class L2:
    @staticmethod
    def handler(samplePPIr, nodeX, nodeY, scoringMethod):
        if scoringMethod == "commonNeighbor":
            score = len(samplePPIr[nodeX]&samplePPIr[nodeY])
        return score

class CAR:
    @staticmethod
    def CRA(samplePPIr, nodeX, nodeY):
        score, cn = 0, samplePPIr[nodeX]&samplePPIr[nodeY]
        for node in cn: score += len(samplePPIr[node]&cn)/len(samplePPIr[node])
        return score

    @staticmethod
    def CAR(samplePPIr, nodeX, nodeY):
        score, cn = 0, samplePPIr[nodeX]&samplePPIr[nodeY]
        for node in cn: score += len(samplePPIr[node]&cn)/2
        return len(cn)*score

    @staticmethod
    def CH2_L3(samplePPIr, nodeX, nodeY):
        uvPair, _, _ = L3.get_uv(nodeX, nodeY, samplePPIr, uvJoin=True)
        U, V = set([uv[0] for uv in uvPair]), set([uv[1] for uv in uvPair])
        localCommunity = U|V
        score = 0
        for [u, v] in uvPair:
            numerator = math.sqrt((1+len(samplePPIr[u]&localCommunity))*(1+len(samplePPIr[v]&localCommunity)))
            denominator = math.sqrt((1+len(samplePPIr[u]-localCommunity-{nodeX, nodeY}))*(1+len(samplePPIr[v]-localCommunity-{nodeX, nodeY})))
            score += numerator/denominator
        return score

    @staticmethod
    def handler(samplePPIr, nodeX, nodeY, scoringMethod):
        if scoringMethod == 'CRA':
            score = CAR.CRA(samplePPIr, nodeX, nodeY)
        if scoringMethod == 'CAR':
            score = CAR.CAR(samplePPIr, nodeX, nodeY)
        if scoringMethod == "CH2_L3":
            score = CAR.CH2_L3(samplePPIr, nodeX, nodeY)
        return score

class L3E:
    jaccard_index = lambda A, B: len(A&B) / len(A|B)
    jaccard_index_alt = lambda A, B: len(A&B) / (len(A|B)-1) if (len(A|B)-1) != 0 else 0
    CRA_index = lambda A, B: len(A&B)/len(A)
    CRA_index_alt = lambda A, B: len(A&B)/(len(A)-1) if (len(A)-1) != 0 else 0
    outer_index = lambda a, B, N: len(B)/len(N[a]) # same for all scoringMethod (except simpCount) because $B \in N(a)$
    func_map = {"f2": jaccard_index, "f1": CRA_index, "f1Alt": CRA_index_alt, "f2Alt": jaccard_index_alt}
    
    @staticmethod
    def handler(samplePPIr, nodePairs, scoringMethod, testmode, logging):
        # scoringMethod structure: L3E_{function name}
        sim_index = None
        for func in L3E.func_map:
            if func == scoringMethod.split("_")[1]: sim_index = L3E.func_map[func]
        # note that for CRA_index, the input order of the arg to sim_index matters
        innerFunc = lambda u,v,U,V,x,y,N : sim_index(N[v],N[x]) * sim_index(N[u],N[y]) * sim_index(N[u],V) * sim_index(N[v],U)
        outerFunc = lambda U,V,x,y,N : L3E.outer_index(x,U,N) * L3E.outer_index(y,V,N)

        scores, predictedPPIbrs = [], []
        count, lastCount, total, avgTime = 0, 0, len(nodePairs), 0
        for nodePair in nodePairs:
            startTime = time.time()
            if nodePair[1] in samplePPIr[nodePair[0]] and not testmode: continue
            nodeX, nodeY = nodePair[0], nodePair[1]
            uvPair, candidateUs, candidateVs = L3.get_uv(nodeX, nodeY, samplePPIr, uvJoin=True)
            nodeUs, nodeVs = set([uv[0] for uv in uvPair]), set([uv[1] for uv in uvPair])
            score = 0
            for [u, v] in uvPair:
                score += innerFunc(u,v,nodeUs,nodeVs,nodeX,nodeY,samplePPIr)
            score *= outerFunc(nodeUs,nodeVs,nodeX,nodeY,samplePPIr)
            scores.append(score)
            predictedPPIbrs.append(nodePair)
            if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)
        return scores, predictedPPIbrs
        

def _PPILinkPred(nodePairs, samplePPIr, scoringMethod, logging=False, testmode=False):
    scores, predictedPPIbrs = [], []
    count, lastCount, total, avgTime = 0, 0, len(nodePairs), 0
    startTime = time.time()

    if scoringMethod in ns.L3Based:
        for nodePair in nodePairs:
            startTime = time.time()
            if nodePair[1] in samplePPIr[nodePair[0]] and not testmode: continue
            nodeX, nodeY = nodePair[0], nodePair[1]
            score = L3.handler(samplePPIr, nodeX, nodeY, scoringMethod)
            scores.append(score)
            predictedPPIbrs.append(nodePair)
            if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)

    elif scoringMethod in ns.L2Based:
        for nodePair in nodePairs:
            startTime = time.time()
            if nodePair[1] in samplePPIr[nodePair[0]] and not testmode: continue
            nodeX, nodeY = nodePair[0], nodePair[1]
            score = L2.handler(samplePPIr, nodeX, nodeY, scoringMethod)
            scores.append(score)
            predictedPPIbrs.append(nodePair)
            if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)

    elif scoringMethod in ns.CARBased:
        for nodePair in nodePairs:
            startTime = time.time()
            if nodePair[1] in samplePPIr[nodePair[0]] and not testmode: continue
            nodeX, nodeY = nodePair[0], nodePair[1]
            score = CAR.handler(samplePPIr, nodeX, nodeY, scoringMethod)
            scores.append(score)
            predictedPPIbrs.append(nodePair)
            if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)

    elif "L3E" in scoringMethod:
        scores, predictedPPIbrs = L3E.handler(samplePPIr, nodePairs, scoringMethod, testmode, logging)

    elif scoringMethod in ns.controlBased:
        for nodePair in nodePairs:
            startTime = time.time()
            if nodePair[1] in samplePPIr[nodePair[0]] and not testmode: continue
            nodeX, nodeY = nodePair[0], nodePair[1]
            score = control.handler(samplePPIr, nodeX, nodeY, scoringMethod)
            scores.append(score)
            predictedPPIbrs.append(nodePair)
            if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)

    # print("core finished")
    return scores, predictedPPIbrs

def _multiCore_handler(args, iterable):
    (nodePairs, splitStartIndex, splitEndIndex, samplePPIr, scoringMethod, logging, testmode, PPIresQ) = args
    nodePairs = nodePairs[splitStartIndex[iterable]:splitEndIndex[iterable]]
    logging = logging[iterable]
    scores, predictedPPIbrs = _PPILinkPred(nodePairs, samplePPIr, scoringMethod, logging, testmode)
    PPIresQ.put([predictedPPIbrs, scores])
    return

def multiCore_PPILinkPred(samplePPIbr, scoringMethod, coreNo, topNo=None, logging=False, testmode=False, nodePairs=None):

    samplePPIr = ns.BRToRelat(ns.toDualBR(samplePPIbr), rSet=True)
    sampleNodes = ns.BRToNode(samplePPIbr)
    if nodePairs is None: nodePairs = list(combinations(sampleNodes, 2))

    if scoringMethod == "random":
        candidatePPIs = ns.pStr_arr(set(ns.arr_pStr(nodePairs))-set(ns.arr_pStr(ns.toDualBR(samplePPIbr))))
        rn.shuffle(candidatePPIs)
        if topNo is None: topNo = len(candidatePPIs)
        return candidatePPIs[0:topNo], [1 for i in range(topNo)]

    splitStartIndex = [i*math.floor(len(nodePairs)/coreNo) for i in range(0, coreNo)] # both splitting is correct
    splitEndIndex = [(i+1)*math.floor(len(nodePairs)/coreNo) if i != coreNo-1 else len(nodePairs) for i in range(0, coreNo)]
    mgr = Manager()
    PPIresQ = mgr.Queue()
    if logging: logging = [True if i == 0 else False for i in range(coreNo)]
    else: logging = [False for i in range(coreNo)]
    args = (nodePairs, splitStartIndex, splitEndIndex, samplePPIr, scoringMethod, logging, testmode, PPIresQ)
    func = partial(_multiCore_handler, args)
    with Pool(coreNo) as p:
        p.map(func, [i for i in range(coreNo)])
    if logging: print("\n")
    mergedScores, mergedPPIbrs = [], []
    PPIresL = [PPIresQ.get() for i in range(coreNo)]
    for [predictedPPIbr, scores] in PPIresL:
        mergedScores += scores
        mergedPPIbrs += predictedPPIbr

    sortedPPIbrs, sortedScores = hr.sort_key_val(mergedPPIbrs, mergedScores)
    if topNo is None: topNo = len(sortedPPIbrs)
    topPredPPIbrs = sortedPPIbrs[0:topNo]
    topScores = sortedScores[0:topNo]
    return topPredPPIbrs, topScores
    
def get_sliding_prec(fullPPIbr, topPredPPIbrs, loopRange, logging):
    minI, maxI = loopRange[0], loopRange[1]
    fullPPIbr = set(ns.arr_pStr(ns.toDualBR(fullPPIbr)))
    precTop, precBot = len(set(ns.arr_pStr(topPredPPIbrs[0:minI]))&fullPPIbr), minI
    precTops = []
    precs = [precTop/precBot]
    count, lastCount, total, avgTime, startTime = minI+1, 0, maxI, 0, time.time()
    for i in range(minI+1, maxI):
        if ns.br_str(topPredPPIbrs[i]) in fullPPIbr or ns.br_str(topPredPPIbrs[i][::-1]) in fullPPIbr: precTop += 1
        precTops.append(precTop)
        if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)
    precTops = np.asarray(precTops)
    precBots = np.asarray([i for i in range(minI+1, maxI)])
    precs += list(np.divide(precTops, precBots))
    return precs

def get_sliding_rec(fullPPIbr, samplePPIbrs, topPredPPIbrs, loopRange, logging):
    minI, maxI = loopRange[0], loopRange[1]
    relevant = set(ns.arr_pStr(fullPPIbr))-set(ns.arr_pStr(ns.toDualBR(samplePPIbrs)))
    fullPPIbr = set(ns.arr_pStr(ns.toDualBR(fullPPIbr)))
    recTop, recBot = len(set(ns.arr_pStr(topPredPPIbrs[0:minI]))&fullPPIbr), len(relevant)
    recTops = []
    recs = [recTop/recBot]
    count, lastCount, total, avgTime, startTime = minI+1, 0, maxI, 0, time.time()
    for i in range(minI+1, maxI):
        if ns.br_str(topPredPPIbrs[i]) in fullPPIbr or ns.br_str(topPredPPIbrs[i][::-1]) in fullPPIbr: recTop += 1
        recTops.append(recTop)
        if logging: count, lastCount, avgTime = helperFunc.logging(count, lastCount, total, avgTime, startTime)
    recTops = np.asarray(recTops)
    recBots = np.asarray([recBot for i in range(minI+1, maxI)])
    recs += list(np.divide(recTops, recBots))
    return recs

def get_sliding_fPosR(fullPPIbr, samplePPIbrs, topPredPPIbrs, loopRange, logging):
    minI, maxI = loopRange[0], loopRange[1]
    fullPPIbr = set(ns.arr_pStr(ns.toDualBR(fullPPIbr)))
    sampleNodes = ns.BRToNode(samplePPIbrs)
    realNeg = len(list(combinations(sampleNodes, 2)))-len(fullPPIbr)
    fPosNum = len(set(ns.arr_pStr(topPredPPIbrs[0:minI]))-fullPPIbr)
    fPoss = []
    fPosRs = [fPosNum/realNeg]
    for i in range(minI+1, maxI):
        if ns.br_str(topPredPPIbrs[i]) not in fullPPIbr or ns.br_str(topPredPPIbrs[i][::-1]) not in fullPPIbr: fPosNum += 1
        fPoss.append(fPosNum)
    fPoss = np.asarray(fPoss)
    realNegs = np.asarray([realNeg for i in range(minI+1, maxI)])
    fPosRs += list(np.divide(fPoss, realNegs))
    return fPosRs

def rocMap_handler(args, iterable):
    (perTagNo, tags, predPPIbr, samplePPIbr, fullPPIbr, logging, resQ) = args
    logging = logging[iterable]
    partialROCMap = {}
    for tagNo in perTagNo[iterable]:
        loopRange = (1, len(predPPIbr[tagNo]))
        partialROCMap[tags[tagNo]] = {}
        partialROCMap[tags[tagNo]]["rocX"] = get_sliding_fPosR(fullPPIbr[tagNo], samplePPIbr[tagNo], predPPIbr[tagNo], loopRange, logging)
        partialROCMap[tags[tagNo]]["rocY"] = get_sliding_rec(fullPPIbr[tagNo], samplePPIbr[tagNo], predPPIbr[tagNo], loopRange, logging)
    resQ.put(partialROCMap)
    return

def rocMap_multiCore(tags, predPPIbr, samplePPIbr, fullPPIbr, coreNo, logging=False):
    tagNo, perTagNo = [i for i in range(len(tags))], [[] for i in range(coreNo)]
    for i in range(math.ceil(len(tags)/coreNo)):
        for j in range(coreNo*i, min(coreNo*(i+1), coreNo*i+(len(tags)-coreNo*i))):
            perTagNo[j-coreNo*i].append(tagNo[j])
    mgr = Manager()
    resQ = mgr.Queue()
    if logging: logging = [True if i == 0 else False for i in range(coreNo)]
    else: logging = [False for i in range(coreNo)]
    args = (perTagNo, tags, predPPIbr, samplePPIbr, fullPPIbr, logging, resQ)
    func = partial(rocMap_handler, args)
    with Pool(coreNo) as p:
        p.map(func, [i for i in range(coreNo)])
    rocMap = {}
    for i in range(coreNo): rocMap.update(resQ.get())
    return rocMap

def precRecMap_handler(args, iterable):
    (perTagNo, tags, predPPIbr, samplePPIbr, fullPPIbr, logging, resQ) = args
    logging = logging[iterable]
    partialPrecRecMap = {}
    for tagNo in perTagNo[iterable]:
        loopRange = (1, len(predPPIbr[tagNo]))
        partialPrecRecMap[tags[tagNo]] = {}
        partialPrecRecMap[tags[tagNo]]["prec"] = get_sliding_prec(fullPPIbr[tagNo], predPPIbr[tagNo], loopRange, logging)
        partialPrecRecMap[tags[tagNo]]["rec"] = get_sliding_rec(fullPPIbr[tagNo], samplePPIbr[tagNo], predPPIbr[tagNo], loopRange, logging)
    resQ.put(partialPrecRecMap)
    return

def precRecMap_multiCore(tags, predPPIbr, samplePPIbr, fullPPIbr, coreNo, logging=False):
    tagNo, perTagNo = [i for i in range(len(tags))], [[] for i in range(coreNo)]
    for i in range(math.ceil(len(tags)/coreNo)):
        for j in range(coreNo*i, min(coreNo*(i+1), coreNo*i+(len(tags)-coreNo*i))):
            perTagNo[j-coreNo*i].append(tagNo[j])
    mgr = Manager()
    resQ = mgr.Queue()
    if logging: logging = [True if i == 0 else False for i in range(coreNo)]
    else: logging = [False for i in range(coreNo)]
    args = (perTagNo, tags, predPPIbr, samplePPIbr, fullPPIbr, logging, resQ)
    func = partial(precRecMap_handler, args)
    with Pool(coreNo) as p:
        p.map(func, [i for i in range(coreNo)])
    precRecMap = {}
    for i in range(coreNo): precRecMap.update(resQ.get())
    return precRecMap

if __name__ == "__main__":
    pass
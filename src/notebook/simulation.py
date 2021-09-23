import sys, os
sys.path.append('..')
import PPILinkPred as ppiLPred
import traversalHelper as tr
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score
import numpy as np
from collections import defaultdict
import math
import random as rand
from statistics import median
from statistics import mean
import multiprocessing
from multiprocessing import Pool, Manager
from functools import partial
from itertools import combinations

class old:
    # deterministics 
    def L3Bridge_graph(nPath, nBridge=None):
        # bridging v to u
        if nBridge is None: nBridge = (nPath**2)-nPath # default bridge all possible paths
        # example, given 3 paths, there are 3*3 edges, 3 edges already exist due to P4, so n*n-n
        if nBridge > (nPath**2)-nPath: raise ValueError("invalid number of bridges")
        
        N = simpleL3_graph(nPath)
        
        edgeCount = 0
        for i in range(nPath): # bridge u to other v
            curNode = 'u'+str(i)
            for j in range(nPath):
                if edgeCount == nBridge: break
                if 'v'+str(j) in N[curNode]: continue
                N[curNode].add('v'+str(j))
                N['v'+str(j)].add(curNode)
                edgeCount += 1
            if edgeCount == nBridge: break

        return N

    # deterministics, missing some incompatible edges
    def L3cluster_graph(nPath, nBridge=None):
        # bridging u to y or v to x
        # number of path * 2 >= nBridge
        # if nBridge not provided, flood the cluster
        if nBridge is None: nBridge = nPath*2
        if nBridge > nPath*2: raise ValueError("invalid number of bridges") 
        
        N = simpleL3_graph(nPath)
        
        for i in range(nPath*2):
            if i >= nBridge: break
                
            if i % 2 == 0:
                node = 'u'+str(int(i/2)) # connect u node to y node to "bridge" a cluster
                N[node].add('y')
                N['y'].add(node)
            else:
                node = 'v'+str(int(i/2)) # connect v node to x node to "bridge" a cluster
                N[node].add('x')
                N['x'].add(node)

        return N



class init_graph:
    @staticmethod
    def cn_graph(graphSize):
        nodeBs = ['a'+str(i) for i in range(graphSize)]
        br = []
        for node in nodeBs: br += [['x', node], [node, 'x'], ['y', node], [node, 'y']]
        N = tr.Helper.binary_to_relation(br, rSet=True) # neighbor func as a set 
        return N # neighbor func represents graph

    @staticmethod
    def ideal_L3_graph(graphSize):
        br = []
        uvs = [("u"+str(i), "v"+str(j)) for i in range(graphSize) for j in range(graphSize)]
        for i in range(graphSize): br.append(('x', 'u'+str(i)))
        for i in range(graphSize): br.append(('y', 'v'+str(i))) 
        br += uvs
        N = tr.Helper.binary_to_relation(br, rSet=True)
        return N

    @staticmethod
    def simple_L3_graph(graphSize):
        '''
        generate a graph consists of n simple L3 paths between X Y
        '''
        nodeUs = ['u'+str(i) for i in range(graphSize)]
        nodeVs = ['v'+str(i) for i in range(graphSize)]
        br = []
        for node in nodeUs: br += [['x', node], [node, 'x']]
        for node in nodeVs: br += [['y', node], [node, 'y']]
        for i in range(graphSize): br += [[nodeUs[i], nodeVs[i]], [nodeVs[i], nodeUs[i]]]
        N = tr.Helper.binary_to_relation(br, rSet=True) # neighbor func as a set
        return N # neighbor func represents graph


class mod_graph:
    '''
    @graphSize: number of u = number of v = number of x edge = number of y edge
    @n: number of edges to add / remove
    '''
    max_edge_num = {
    'removeComp_idealL3': lambda graphSize : graphSize**2-graphSize,
    'addComp_simpleL3': lambda graphSize : graphSize**2-graphSize,
    'addIncompA': lambda graphSize : graphSize*(graphSize-1)+graphSize*2,
    'addIncompB': lambda graphSize : graphSize*(graphSize-1),
    'addIncompC': lambda graphSize : graphSize*2
    }

    @staticmethod
    def removeComp_idealL3(graphSize, n, init_graph_func):
        '''
        removing compatible edges between u nodes & v nodes
        must keep edges such that any u (v) is adjacent to at least one v (u)
        '''
        if n > mod_graph.max_edge_num['removeComp_idealL3'](graphSize):
            raise ValueError("invalid number of compatible edges")

        N = init_graph_func(graphSize) # has to be idealL3
        candidates = [("u"+str(i), "v"+str(j)) for i in range(graphSize) for j in range(graphSize) if i != j]
        for candidate in rand.sample(candidates, n):
            N[candidate[0]].remove(candidate[1])
            N[candidate[1]].remove(candidate[0])
        return N

    @staticmethod
    def addComp_simpleL3(graphSize, n, init_graph_func):
        '''
        add edges between u and v
        '''
        if n > mod_graph.max_edge_num['addComp_simpleL3'](graphSize):
            raise ValueError("invalid number of compatible edges")

        N = init_graph_func(graphSize) # has to be simpleL3
        # generate a list of possible paths
        candidates = [("u"+str(i), "v"+str(j)) for i in range(graphSize) for j in range(graphSize) if i != j]
        for candidate in rand.sample(candidates, n):
            N[candidate[0]].add(candidate[1])
            N[candidate[1]].add(candidate[0])
        return N

    @staticmethod
    def addIncompA(graphSize, n, init_graph_func):
        '''
        add incompatible edges between pair of u or pair of v or x v or y u
        '''
        if n > mod_graph.max_edge_num['addIncompA'](graphSize):
            raise ValueError("invalid number of incompatible edges")

        N = init_graph_func(graphSize)
        candidates = [("u"+str(i), "u"+str(j)) for (i,j) in list(combinations([i for i in range(graphSize)], 2))]
        candidates += [("v"+str(i), "v"+str(j)) for (i,j) in list(combinations([i for i in range(graphSize)], 2))]
        candidates += [('x', 'v'+str(i)) for i in range(graphSize)] + [('y', 'u'+str(i)) for i in range(graphSize)]
        for candidate in rand.sample(candidates, n):
            N[candidate[0]].add(candidate[1])
            N[candidate[1]].add(candidate[0])
        return N

    @staticmethod
    def addIncompB(graphSize, n, init_graph_func):
        '''
        add incompatible edges between pair of u or pair of v
        '''
        if n > mod_graph.max_edge_num['addIncompB'](graphSize):
            raise ValueError("invalid number of incompatible edges")

        N = init_graph_func(graphSize)
        candidates = [("u"+str(i), "u"+str(j)) for (i,j) in list(combinations([i for i in range(graphSize)], 2))]
        candidates += [("v"+str(i), "v"+str(j)) for (i,j) in list(combinations([i for i in range(graphSize)], 2))]
        for candidate in rand.sample(candidates, n):
            N[candidate[0]].add(candidate[1])
            N[candidate[1]].add(candidate[0])
        return N

    @staticmethod
    def addIncompC(graphSize, n, init_graph_func):
        '''
        add incompatible edges between x v or y u
        '''
        if n > mod_graph.max_edge_num['addIncompC'](graphSize):
            raise ValueError("invalid number of incompatible edges")

        N = init_graph_func(graphSize)
        candidates = [('x', 'v'+str(i)) for i in range(graphSize)] + [('y', 'u'+str(i)) for i in range(graphSize)]
        for candidate in rand.sample(candidates, n):
            N[candidate[0]].add(candidate[1])
            N[candidate[1]].add(candidate[0])
        return N


# placing exactly here is important
class ns:
    init_graph_map = {
        "cn_graph": init_graph.cn_graph,
        "ideal_L3_graph": init_graph.ideal_L3_graph,
        "simple_L3_graph": init_graph.simple_L3_graph
    }

    mod_graph_map = {
        "removeComp_idealL3": mod_graph.removeComp_idealL3,
        "addComp_simpleL3": mod_graph.addComp_simpleL3,
        "addIncompA": mod_graph.addIncompA,
        "addIncompB": mod_graph.addIncompB,
        "addIncompC": mod_graph.addIncompC
        }


def sim_core(args, iterable):
    (scores, graphSize, n, method, graph_func, simulating, coreNo, trials) = args
    scores_loc = [{} for i in range(trials)]
    startInd, endInd = int(n/coreNo*iterable), int(n/coreNo*(iterable+1))
    for trial in range(trials):
        for n in range(startInd, endInd):
            N = ns.mod_graph_map[simulating](
                graphSize, n, ns.init_graph_map[graph_func])
            score, _ = ppiLPred._PPILinkPred([['x', 'y']], N, scoringMethod=method)
            scores_loc[trial][n] = score[0]
    scores.put(scores_loc)
    return

def simulation(method, graphSize, saturateRatio, graph_func, simulating, coreNo=10, trials=10, normalizing=True):
    '''
    @graphSize: see class mod_graph
    @saturateRatio: ratio of edges number to add w.r.t to the valid maximum
    @graph_func: (STRING) func name, see func_map in class init_graph
    @simulating: (STRING) func name, see func_map in class mod_graph
    '''
    if (simulating == 'removeComp_idealL3' and graph_func != 'ideal_L3_graph') or (
        simulating == 'addComp_simpleL3' and graph_func != "simple_L3_graph"):
        raise ValueError("invalid simulating function for the init graph type")

    maxEdgeNum = mod_graph.max_edge_num[simulating](graphSize)
    n = int(maxEdgeNum*saturateRatio)+1 # off-by-one to offset end of for-loop

    if coreNo > 1:
        mgr = Manager()
        scores = mgr.Queue()
        args = (scores, graphSize, n, method, graph_func, simulating, coreNo, trials)
        func = partial(sim_core, args)
        with Pool(coreNo) as p:
            p.map(func, [i for i in range(coreNo)])
        res = [scores.get() for i in range(coreNo)]
        mergedRes = [{} for i in range(trials)]

        for scoresDict in res:
            for trial in range(trials):
                mergedRes[trial].update(scoresDict[trial])
        scores = [np.asarray([res[n] for n in sorted(res)]) for res in mergedRes]

        if normalizing:
            for i in range(len(scores)):
                scores[i] = (scores[i]-np.min(scores[i]))/(np.max(scores[i])-np.min(scores[i]))
            scores[i][scores[i] == np.inf] = 0
    else:
        scores = []
        for trial in range(trials):
            scores.append([])
            for j in range(n):
                N = ns.mod_graph_map[simulating](
                    graphSize, j, ns.init_graph_map[graph_func])
                score, _ = ppiLPred._PPILinkPred([['x', 'y']], N, scoringMethod=method)
                scores[-1].append(score[0])

            if normalizing:
                scores[-1] = np.asarray(scores[-1])
                scores[-1] = (scores[-1]-np.min(scores[-1]))/(np.max(scores[-1])-np.min(scores[-1]))
                scores[-1][scores[-1] == np.inf] = 0

    return scores

if __name__ == "__main__":
    scoresA = simulation(method="countP4", graphSize=10, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="removeComp_idealL3", coreNo=1, trials=2)
    scoresB = simulation(method="countP4", graphSize=10, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="removeComp_idealL3", coreNo=2, trials=2)
    # avoid float precision inaccuracy
    scoresA = [np.around(i, 4) for i in np.asarray(scoresA)[0]]
    scoresB = [np.around(i, 4) for i in np.asarray(scoresB)[0]]
    assert scoresA[0] == scoresB[0]
    assert scoresA[-1] == scoresB[-1]    
    print(scoresA)

    quit()
    # unit test
    # only assert start score and end score since the rest are randomized
    scoresA = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="removeComp_idealL3", coreNo=1, trials=2, normalizing=False)
    scoresB = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="removeComp_idealL3", coreNo=2, trials=2, normalizing=False)
    # avoid float precision inaccuracy
    scoresA = [np.around(i, 4) for i in np.asarray(scoresA)[0]]
    scoresB = [np.around(i, 4) for i in np.asarray(scoresB)[0]]
    assert scoresA[0] == scoresB[0]
    assert scoresA[-1] == scoresB[-1]
    
    scoresA = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="addIncompA", coreNo=1, trials=1, normalizing=False)
    scoresB = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="addIncompA", coreNo=2, trials=1, normalizing=False)
    # avoid float precision inaccuracy
    scoresA = [np.around(i, 4) for i in np.asarray(scoresA)[0]]
    scoresB = [np.around(i, 4) for i in np.asarray(scoresB)[0]]
    assert scoresA[0] == scoresB[0]
    assert scoresA[-1] == scoresB[-1]

    scoresA = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="addIncompB", coreNo=1, trials=2)
    scoresB = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="addIncompB", coreNo=2, trials=2)
    # avoid float precision inaccuracy
    scoresA = [np.around(i, 4) for i in np.asarray(scoresA)[0]]
    scoresB = [np.around(i, 4) for i in np.asarray(scoresB)[0]]
    assert scoresA[0] == scoresB[0]
    assert scoresA[-1] == scoresB[-1]

    scoresA = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="addIncompC", coreNo=1, trials=2)
    scoresB = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="ideal_L3_graph", simulating="addIncompC", coreNo=2, trials=2)
    # avoid float precision inaccuracy
    scoresA = [np.around(i, 4) for i in np.asarray(scoresA)[0]]
    scoresB = [np.around(i, 4) for i in np.asarray(scoresB)[0]]
    assert scoresA[0] == scoresB[0]
    assert scoresA[-1] == scoresB[-1]

    scoresA = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="simple_L3_graph", simulating="addComp_simpleL3", coreNo=1, trials=2)
    scoresB = simulation(method="L3E1_f1", graphSize=3, saturateRatio=1,
    graph_func="simple_L3_graph", simulating="addComp_simpleL3", coreNo=2, trials=2)
    # avoid float precision inaccuracy
    scoresA = [np.around(i, 4) for i in np.asarray(scoresA)[0]]
    scoresB = [np.around(i, 4) for i in np.asarray(scoresB)[0]]
    assert scoresA[0] == scoresB[0]
    assert scoresA[-1] == scoresB[-1]

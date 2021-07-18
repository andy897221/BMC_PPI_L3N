import math
from itertools import combinations
# see unittest.py to see the verification of L3E

class helper:
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


class L3:
    @staticmethod
    def L3(N):
        scores, predictedPPIs = [], []
        allNodePairs = list(combinations(list(N.keys()), 2))
        for nodePair in allNodePairs:
            if nodePair[1] in N[nodePair[0]]: continue
            x, y = nodePair[0], nodePair[1]
            uvPair, _, _ = helper.get_uv(x, y, N)
            score = 0
            for [u, v] in uvPair:
                if math.sqrt(len(N[u])*len(N[v])) == 0: continue
                score += 1/math.sqrt(len(N[u])*len(N[v]))
            scores.append(score)
            predictedPPIs.append(nodePair)
        return scores, predictedPPIs

class L3E:
    sim_f1 = lambda A, B: len(A&B)/len(A) # similarity function: simple ratio
    sim_f2 = lambda A, B: len(A&B)/len(A|B) # similarity function: jaccard index

    @staticmethod
    def L3E(N, sim_name):
        if sim_name == 'f1': sim_f = L3E.sim_f1
        elif sim_name == 'f2': sim_f = L3E.sim_f2
        else: raise AttributeError("similarity function {} does not exist".format(sim_name))
        scores, predictedPPIs = [], []
        allNodePairs = list(combinations(list(N.keys()), 2))
        for nodePair in allNodePairs:
            if nodePair[1] in N[nodePair[0]]: continue
            x, y = nodePair[0], nodePair[1]
            uvPair, candidateUs, candidateVs = helper.get_uv(x, y, N, uvJoin=True)
            U, V = set([uv[0] for uv in uvPair]), set([uv[1] for uv in uvPair])
            score = 0
            for [u, v] in uvPair:
                score += sim_f(N[v],N[x]) * sim_f(N[u],N[y]) * sim_f(N[u],V) * sim_f(N[v],U)
            score *= sim_f(N[x],U) * sim_f(N[y],V)
            scores.append(score)
            predictedPPIs.append(nodePair)
        return scores, predictedPPIs
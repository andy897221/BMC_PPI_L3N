import sys
import pandas as pd
import numpy as np
import core
from collections import defaultdict

def binary_to_neighborhood(br):
    relation = defaultdict(set)
    for nodeArr in br:
        relation[nodeArr[0]].add(nodeArr[1])
        relation[nodeArr[1]].add(nodeArr[0])
    return relation

if __name__ == "__main__":
    ppiFile, savePath, method = sys.argv[1], sys.argv[2], sys.argv[3]

    ppi_df = pd.read_csv(ppiFile, sep="\t", header=None)
    ppi_br = ppi_df[[0,1]].values.tolist()
    N = binary_to_neighborhood(ppi_br)

    scoreArgs = {}
    if method == "L3Ef1":
        sortedScores, sortedPPIs = core.L3E.L3E(N, "f1")
    elif method == "L3Ef2":
        sortedScores, sortedPPIs = core.L3E.L3E(N, "f2")
    elif method == "L3":
        sortedScores, sortedPPIs = core.L3.L3(N)
    else:
        print("link predictor '{}' is not supported / invalid.".format(method))
        exit()

    sortedPPIs = np.asarray(sortedPPIs).transpose()
    save_df = pd.DataFrame({"proteinA": sortedPPIs[0], "proteinB": sortedPPIs[1], "score": sortedScores})
    save_df.to_csv(savePath, sep="\t", index=False)
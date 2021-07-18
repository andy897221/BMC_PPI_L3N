import numpy as np
from operator import itemgetter
import sys

def get_item_I(arr, item):
    # unique item
    return np.where(np.asarray(arr) == item)[0][0]


def key_sorted_by_val(d, lowToHigh=False):
    if lowToHigh:
        sortedD = sorted(d.items(), key=itemgetter(1))
    else:
        sortedD = sorted(d.items(), key=itemgetter(1), reverse=True)
    sortedKey = [k for k, _ in sortedD]
    return sortedKey

def sort_key_val(keys, vals, lowToHigh=False):
    valsMap = {}
    for i in range(0, len(vals)): valsMap[i] = vals[i]
    sortedMap = key_sorted_by_val(valsMap)
    sortedKeys = [keys[i] for i in sortedMap]
    sortedVals = [vals[i] for i in sortedMap]
    return sortedKeys, sortedVals

if __name__ == "__main__":
    pass
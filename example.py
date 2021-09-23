import sys
sys.path.append("./src")
import core

# consider such a PPI network
#  A
# / \
#B  E
#|  |
#C  F
#\ /
# D
# expressed as a neighborhood dictionary
PPIr = {"A": {"B", "E"}, "B": {"A", "C"}, "C": {"B", "D"}
, "E": {"A", "F"}, "F": {"E", "D"}, "D": {"C", "F"}}

# to generate all candidate edges with scores using the original L3 link predictor
sortedScores, sortedPPIs = core.L3.L3(PPIr)
print("L3:\n {} \n {}\n".format(sortedPPIs, sortedScores))

# to generate all candidate edges with scores using the L3E(f1) link predictor
sortedScores, sortedPPIs = core.L3E.L3E(PPIr, "f1")
print("L3EPrime(f1):\n {} \n {}\n".format(sortedPPIs, sortedScores))

# to generate all candidate edges with scores using the L3E(f1) link predictor
sortedScores, sortedPPIs = core.L3E.L3E(PPIr, "f1Alt")
print("L3E(f1):\n {} \n {}\n".format(sortedPPIs, sortedScores))

# in this example PPI network, both rank candidate PPIs the same
# for examples of practical scenarios where L3E performs better, they are documented in our paper
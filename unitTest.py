from itertools import combinations
import core
import sys
import unittest
sys.path.append("./")
sys.path.append("./src")
import PPILinkPred as ppiLPred

class testData:
    '''
    consider such a PPI network, a example CN structure (shared 3 nodes, one dangling node)
          A
        / | \
        B C D
        \ |
          E

    consider such a PPI network, a example L3 structure (2 P_{4})
         A
        / \
       B   E
       |   |
       C   F
        \ /
         D
    '''

    @staticmethod
    def getTestA():
        # testA: example L3 structure (above)
        return {"A": {"B", "E"}, "B": {"A", "C"}, "C": {"B", "D"}
        , "E": {"A", "F"}, "F": {"E", "D"}, "D": {"C", "F"}}

    @staticmethod
    def getTestB():
        # testB: testA with additional compatible edge
        testA_PPIr = testData.getTestA()
        testA_PPIr["E"].add("C")
        return testA_PPIr

    @staticmethod
    def getTestC():
        # testC: testC with additional incompatible edge
        testA_PPIr = testData.getTestA()
        testA_PPIr["E"].add("B")
        return testA_PPIr

    @staticmethod
    def getTestD():
        # testD: example CN structure (above)
        return {"A": {"B", "C", "D"}, "B": {"A", "E"}
        , "C": {"A", "E"}, "D": {"A"}, "E": {"B", "C"}}


class unitTestL3E(unittest.TestCase):
    def runTest(self, testcase, allNodePairs):
        self.assertEqual(core.L3E.L3E(testcase, "f1")
        , ppiLPred._PPILinkPred(allNodePairs, testcase, "L3E_f1"))
        self.assertEqual(core.L3E.L3E(testcase, "f2")
        , ppiLPred._PPILinkPred(allNodePairs, testcase, "L3E_f2"))
        self.assertEqual(core.L3E.L3E(testcase, "f1Alt")
        , ppiLPred._PPILinkPred(allNodePairs, testcase, "L3E_f1Alt"))
        self.assertEqual(core.L3E.L3E(testcase, "f2Alt")
        , ppiLPred._PPILinkPred(allNodePairs, testcase, "L3E_f2Alt"))

    def testA(self):
        testcaseA = testData.getTestA()
        allNodePairs = list(combinations(list(testcaseA.keys()), 2))
        self.runTest(testcaseA, allNodePairs)

    def testB(self):
        testcaseB = testData.getTestB()
        allNodePairs = list(combinations(list(testcaseB.keys()), 2))
        self.runTest(testcaseB, allNodePairs)

    def testC(self):
        testcaseC = testData.getTestC()
        allNodePairs = list(combinations(list(testcaseC.keys()), 2))
        self.runTest(testcaseC, allNodePairs)

    def testD(self):
        testcaseD = testData.getTestD()
        allNodePairs = list(combinations(list(testcaseD.keys()), 2))
        self.runTest(testcaseD, allNodePairs)

if __name__ == "__main__":
    # unit test: assert with the L3E used the experiments "./src/PPILinkPred.py" (non-cleaned version)
    unittest.main()
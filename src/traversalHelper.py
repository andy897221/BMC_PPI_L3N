import numpy as np
import sys
from collections import defaultdict
import helper as hr

class Helper:
    def pathNotTraveled(paths, path):
        for i in paths:
            if "".join([str(j) for j in path]) in "".join([str(j) for j in i]): return False
        return True

    def isNotCycle(path, node):
        return not node in path

    def edgeStr(nodeArr):
        # OBSELETE, BECAUSE SOME GENE HAS A SYMBOL OF "-"
        return str(nodeArr[0])+"-"+str(nodeArr[1])

    def item_in_set(i, S):
        for j in i:
            if j in S: return True
        return False

    def is_subset(subS, S):
        count = 0
        for i in subS:
            if i in S: count += 1
        if count == len(subS): return True
        return False

    def brsEntry_to_Str(brs):
        brStr = []
        for br in brs:
            brStr.append([str(node) for node in br])
        return brStr

    def binary_relation_to_node(binary_relation):
        return set(list(np.unique(np.asarray(binary_relation).flatten())))

    def pathStr(path):
        # given a list, serialize the list of elements as str, return a str
        return "\t".join([str(i) for i in path])

    def list_to_pathStrs(paths):
        return [Helper.pathStr(path) for path in paths]

    def pathStrs_to_list(pathStr, isInt=False):
        # given a list of str (the str is a serialized list), de-serialize them, return a list of lists
        if isInt:
            paths = [i.split("\t") for i in pathStr]
            intPaths = []
            for path in paths:
                intPaths.append([])
                for i in path: intPaths[-1].append(int(i))
            return intPaths
        else:
            return [i.split("\t") for i in pathStr]

    def br_to_pathStr(br):
        return "{}\t{}".format(br[0], br[1])

    def get_targets(binary_relation):
        source = [i[0] for i in binary_relation]
        target = [i[1] for i in binary_relation]
        leaf = [i for i in target if i not in source]
        return list(set(leaf))

    def get_sources(binary_relation):
        source = [i[0] for i in binary_relation]
        target = [i[1] for i in binary_relation]
        head = [i for i in source if i not in target]
        return list(set(head))

    def to_dual_binary_relation(binary_relation):
        binary_relation_plus = binary_relation.copy()
        for i in binary_relation:
            binary_relation_plus.append([i[1], i[0]])
        binary_relation_plus = Helper.list_to_pathStrs(binary_relation_plus)
        binary_relation_plus = list(set(binary_relation_plus))
        binary_relation_plus = Helper.pathStrs_to_list(binary_relation_plus)
        return binary_relation_plus

    def to_dual_relation(relation):
        relation_plus = relation.copy()
        for key, nodes in relation.items():
            for node in nodes:
                if node not in relation_plus: relation_plus[node] = set([key])
                else: relation_plus[node].add(key)
        return relation_plus

    # this is probably unstable
    # def dual_to_mono(brs):
    #     monoBR = []
    #     for br in brs:
    #         if br[1] < br[0]: monoBR.append([br[1], br[0]])
    #         else: monoBR.append([br[0], br[1]])
    #     monoBR = Helper.pathStrs_to_list(set(Helper.list_to_pathStrs(monoBR)))
    #     return monoBR

    def filter_self_cycle(brs):
        noCycleBR = []
        for br in brs:
            if br[0] == br[1]: continue
            noCycleBR.append([br[0], br[1]])
        return noCycleBR

    def binary_to_relation(binaryRelat, rSet=False):
        nodes = Helper.binary_relation_to_node(binaryRelat)
        # build dict
        relation = {}
        for node in nodes: relation[node] = set()
        for nodeArr in binaryRelat: relation[nodeArr[0]].add(nodeArr[1])
        for nodeArr in binaryRelat: relation[nodeArr[1]].add(nodeArr[0])
        if not rSet:
            for node in relation: relation[node] = list(relation[node])
        return relation

    def relation_to_binary(relation):
        br = []
        for key, arr in relation.items():
            for item in arr: br.append([key, item])
        return br

    def pathsEntry_to_geneName(entry, paths):
        namedPaths = []
        for path in paths:
            namedPaths.append([entry[node]["first_graphics"] for node in path])
        return namedPaths

    def paths_to_binary_relation(paths, dual=False):
        relation = []
        for path in paths:
            for node in range(0, len(path)-1):
                if [path[node], path[node+1]] not in relation:
                    relation.append([path[node], path[node+1]])
                if [path[node+1], path[node]] not in relation and dual:
                    relation.append([path[node+1], path[node]])
        return relation

    def geneName_to_entry(entrys, geneName):
        ids = []
        for node in entrys:
            if entrys[node]["first_graphics"] == geneName: ids += [node]
        return ids
                
    def entry_to_geneName(entrys, entry):
        return entrys[entry]["first_graphics"]

class Traversal:

    # source: https://www.geeksforgeeks.org/cycles-of-length-n-in-an-undirected-and-connected-graph/
    def _cycle_detect(graph, marked, depth, vert, start, cycles, curCycle, vertexNo):
        marked[vert] = True
        curCycle.append(vert)
        if depth == 0:
            marked[vert] = False
            if graph[vert][start] == 1:
                cycles.append(curCycle.copy())
                curCycle.pop()
                return
            else:
                curCycle.pop()
                return

        for i in range(vertexNo):
            if not marked[i] and graph[vert][i] == 1:
                Traversal._cycle_detect(graph, marked, depth-1, i, start, cycles, curCycle, vertexNo)
        curCycle.pop()
        marked[vert] = False
        return

    def cycle_detect(fullBR, depth=4):
        fullBR = Helper.to_dual_binary_relation(fullBR)
        nodes = list(Helper.binary_relation_to_node(fullBR))
        vertexNo = len(nodes)
        marked = [False] * vertexNo
        curCycle, cycles = [], []

        graph = []
        # construct graph
        relation = Helper.binary_to_relation(fullBR, rSet=True)
        for node in nodes:
            graph.append([1 if rNode in relation[node] else 0 for rNode in nodes])

        for i in range(vertexNo-(depth-1)):
            Traversal._cycle_detect(graph, marked, depth-1, i, i, cycles, curCycle, vertexNo)
            marked[i] = True

        # delete repeating cycles
        nonRepeatCycles = []
        for cycle in cycles:
            if [cycle[0]]+cycle[1:len(cycle)][::-1] in nonRepeatCycles: continue
            nonRepeatCycles.append(cycle)

        # node index to node gene name
        namedCycle = []
        for cycle in nonRepeatCycles:
            namedCycle.append([nodes[node] for node in cycle])
        return namedCycle

    def _dfs_cycle_detect(source, target, relation, visited, path, paths, dist):
        visited[source] = True
        path.append(source)
        if source in relation:
            if len(path) >= dist:
                if path[0] in relation[source]: paths.append(path.copy())
                path.pop()
                visited[source]=False
                return

            for i in relation[source]:
                if not visited[i]:
                    Traversal._dfs_cycle_detect(i, target, relation, visited, path, paths, dist)        
        path.pop()
        visited[source]=False

    def dfs_cycle_detect(nodes, relation, source, dist):
        visited = {}
        for node in nodes: visited[node] = False
        path, paths = [], []
        target = source
        Traversal._dfs_cycle_detect(source, target, relation, visited, path, paths, dist)

        nonRepeatCycles = []
        for cycle in paths:
            if [cycle[0]]+cycle[1:len(cycle)][::-1] in nonRepeatCycles: continue
            nonRepeatCycles.append(cycle)
        return nonRepeatCycles

    def _single_dfs_cycle_detect(source, target, relation, visited, path, dist):
        visited[source] = True
        path.append(source)
        if source in relation:  
            if len(path) >= dist:
                if path[0] in relation[source]:
                    return True
                else:
                    path.pop()
                    visited[source]=False
                    return
            
            found = False
            for i in relation[source]:
                if not visited[i]:
                    found = Traversal._single_dfs_cycle_detect(i, target, relation, visited, path, dist)   
                if found: return True
            if found: return True
        path.pop()
        visited[source]=False

    def single_dfs_cycle_detect(nodes, relation, source, dist):
        visited = {}
        for node in nodes: visited[node] = False
        path = []
        target = source
        Traversal._single_dfs_cycle_detect(source, target, relation, visited, path, dist)
        return path

    def _masterSlave_cycle_detect(source, target, relation, visited, path, paths, dist, numOfCycles, masterToggle):
        visited[source] = True
        path.append(source)
        if source in relation:
            if len(path) >= dist:
                if path[0] in relation[source]:
                    paths.append(path.copy())
                    if len(paths) >= numOfCycles: return True
                path.pop()
                visited[source] = False

            found = False
            weightedRelation = {}
            for node in relation[source]:
                if node in relation: weightedRelation[node] = len(relation[node])
            if masterToggle: rankedrNodes = hr.key_sorted_by_val(weightedRelation)
            else: rankedrNodes = hr.key_sorted_by_val(weightedRelation)[::-1]
            masterToggle = not masterToggle
            for i in rankedrNodes:
                if not visited[i]:
                    found = Traversal._masterSlave_cycle_detect(i, target, relation, visited, path, paths, dist, numOfCycles, masterToggle)
                if found: return True
            if found: return True
        path.pop()
        visited[source]=False

    def masterSlave_cycle_detect(nodes, relation, source, dist, numOfCycles, sourceAsMaster):
        # traversal weighted by degree of nodes, find cycles from the source (defined by the numbers)
        visited = {}
        for node in nodes: visited[node] = False
        path, paths = [], []
        target = source
        masterToggle = not sourceAsMaster
        Traversal._masterSlave_cycle_detect(source, target, relation, visited, path, paths, dist, numOfCycles, masterToggle=masterToggle)
        return paths

    def _single_masterSlave_cycle_detect(source, target, relation, visited, path, dist, masterToggle, rankingStrategy, strateArg):
        visited[source] = True
        path.append(source)
        if source in relation:
            if len(path) >= dist:
                if path[0] in relation[source]:
                    return True
                else:
                    path.pop()
                    visited[source] = False
                    return

            found = False
            rankedrNodes = rankingStrategy(source, relation, masterToggle, strateArg)
            masterToggle = not masterToggle
            for i in rankedrNodes:
                if not visited[i]:
                    found = Traversal._single_masterSlave_cycle_detect(i, target, relation, visited, path, dist
                        , masterToggle, rankingStrategy, strateArg)
                if found: return True
            if found: return True
        path.pop()
        visited[source]=False

    def single_masterSlave_cycle_detect(nodes, relation, source, dist, sourceAsMaster,rankingStrategy, strateArg):
        # traversal weighted by degree of nodes, find cycles from the source (defined by the numbers)
        visited = {}
        for node in nodes: visited[node] = False
        path, paths = [], []
        target = source
        masterToggle = not sourceAsMaster
        Traversal._single_masterSlave_cycle_detect(source, target, relation, visited, path, dist
            , masterToggle=masterToggle, rankingStrategy=rankingStrategy, strateArg=strateArg)
        return path

    def _single_cycle_detect(graph, marked, depth, vert, start, curCycle, vertexNo):
        marked[vert] = True
        curCycle.append(vert)
        if depth == 0:
            marked[vert] = False
            if graph[vert][start] == 1: return True
            else: curCycle.pop()
            return False

        found = False
        for i in range(vertexNo):
            if found: return found
            if not marked[i] and graph[vert][i] == 1:
                found = Traversal._single_cycle_detect(graph, marked, depth-1, i, start, curCycle, vertexNo)
        if found: return found
        curCycle.pop()
        marked[vert] = False
        return

    def single_cycle_detect(fullBR, depth=4):
        fullBR = Helper.to_dual_binary_relation(fullBR)
        nodes = np.asarray(list(Helper.binary_relation_to_node(fullBR)))
        vertexNo = len(nodes)
        marked = [False] * vertexNo
        curCycle = []

        graph = []
        # construct graph
        relation = Helper.binary_to_relation(fullBR, rSet=True)
        for node in nodes:
            graph.append([1 if rNode in relation[node] else 0 for rNode in nodes])

        found = False
        for i in range(vertexNo-(depth-1)):
            found = Traversal._single_cycle_detect(graph, marked, depth-1, i, i, curCycle, vertexNo)
            if found: break
            marked[i] = True

        # node index to node gene name
        namedCycle = [nodes[node] for node in curCycle]
        return namedCycle

    def full_dfs(nodes, relation, binary_relation, sources=None, targets=None):
        # since for dual relation, full_dfs is not possible (cannot find sources / target), a sources / targets parameter can be supplied for such scenario
        if sources is None: sources = Helper.get_sources(binary_relation)
        if targets is None: targets = Helper.get_targets(binary_relation)
        paths = []
        for source in sources:
            for target in targets:
                paths += list(Traversal.source_target_dfs(nodes, relation, source, target))
        # paths = Helper.pathStrs_to_list(paths, isInt=isInt)
        return paths

    def source_dfs(nodes, relation, binary_relation, source, targets=None):
        # since for dual relation, source_dfs is not possible (cannot find target), a targets parameter can be supplied for such scenario
        paths = []
        if targets is None: targets = Helper.get_targets(binary_relation)
        for target in targets:
            paths += list(Traversal.source_target_dfs(nodes, relation, source, target))
        # paths = Helper.pathStrs_to_list(paths, isInt=isInt)
        return paths

    # https://www.geeksforgeeks.org/find-paths-given-source-destination/
    # RecursionError: maximum recursion depth exceeded in comparison
    def _source_target_dfs(source, target, relation, visited, path, paths, dist=None):
        visited[source] = True
        path.append(source)
        if source in relation:
            if dist is not None:
                if len(path) >= dist:
                    path.pop()
                    visited[source]=False
                    return

            for i in relation[source]:
                if not visited[i]:
                    Traversal._source_target_dfs(i, target, relation, visited, path, paths, dist)        
        path.pop()
        visited[source]=False

    def source_target_dfs(nodes, relation, source, target, dist=None):
        visited = {}
        for node in nodes: visited[node] = False
        path, paths = [], []
        Traversal._source_target_dfs(source, target, relation, visited, path, paths, dist)
        return paths

    def _node_dfs_search(source, relation, visited, fullBR, depth=None, steps=0):
        visited[source] = True
        if source in relation:
            if depth is not None:
                if steps == depth: return

            for i in relation[source]:
                if not visited[i]:
                    fullBR.append([source, i])
                    Traversal._node_dfs_search(i, relation, visited, fullBR, depth, steps+1)
        return
    
    def node_dfs_search(nodes, relation, source, depth=None):
        visited = {}
        for node in nodes: visited[node] = False
        fullBR = []
        Traversal._node_dfs_search(source, relation, visited, fullBR, depth)
        return fullBR

    def node_bfs_search(nodes, relation, source, depth=None):
        fullBR = []
        visited = {}
        for node in nodes: visited[node] = False
        queue = []
        queue.append(source)
        visited[source] = True
        currentLvNodes = len(queue)
        nextLvNodes = 0
        steps = 0
        while queue:
            s = queue.pop(0)
            currentLvNodes -= 1
            for rNode in relation[s]:
                if not visited[rNode]:
                    queue.append(rNode)
                    visited[rNode] = True
                    fullBR.append([s, rNode])
                    nextLvNodes += 1
            if currentLvNodes == 0:
                currentLvNodes = nextLvNodes
                nextLvNodes = 0
                steps += 1
                if depth is not None:
                    if steps == depth: break
        return fullBR

    def full_bfs(nodes, relation, sources):
        return Traversal._bfs(nodes, relation, source=None, target=None, fullSources=sources)

    def _bfs(nodes, relation, source=None, target=None, fullSources=None):
        visited = {}
        for node in nodes: visited[node] = False
        
        # build queue
        curPtr, queue, parentPtrs, parentPaths = 0, [], {}, []
        if source is not None:
            queue.append(source)
            parentPtrs[source] = 0
            parentPaths.append([source])
        else:
            for source in fullSources:
                queue.append(source)
                parentPaths.append([source])
                parentPtrs[source] = curPtr
                curPtr += 1
        
        visited[source] = True
        finalizedPaths = []
        while queue:
            paths, curPtr, pathVisited = [], 0, [False for i in parentPaths]
            source = queue.pop(0)
            if source not in relation: continue
            for node in relation[source]:
                if not visited[node]:
                    if target is not None:
                        if target == node: return parentPaths[parentPtrs[source]] + [node]
                    paths.append(parentPaths[parentPtrs[source]] + [node])
                    pathVisited[parentPtrs[source]] = True
                    parentPtrs[node] = curPtr
                    curPtr += 1
                    queue.append(node)
                    visited[node] = True
                    # print(source, node, parentPaths, paths)

            # refresh unvisited path's ptrs
            for path in parentPaths:
                if not pathVisited[parentPtrs[path[-1]]]:
                    paths.append(parentPaths[parentPtrs[path[-1]]])
                    parentPtrs[path[-1]] = curPtr
                    curPtr += 1
            parentPaths = paths.copy()
        if target is None: return parentPaths
        else: return None

    def source_bfs(nodes, relation, source):
        return Traversal._bfs(nodes, relation, source=source)

    def source_target_bfs(nodes, relation, source, target):
        return Traversal._bfs(nodes, relation, source=source, target=target)

    def full_st_bfs(nodes, relation, source, target):
        paths = []
        while True:
            tempPath = Traversal._bfs(nodes, relation, source=source, target=target)
            if tempPaths is None: break
            paths.append(tempPath)
            print(paths)
            # cut relation
            relation = []
            for nodeArr in relation:
                if nodeArr == [tempPath[-2], tempPath[-1]]: continue
                relation.append(nodeArr)
        ##### UN-DEBUGGED, CANNOT BE USED YET
        return paths

    def dfsPaths_to_geneName(dfsPaths, entry):
        dfsGeneName = []
        for path in dfsPaths:
            dfsGeneName.append([entry[i]['first_graphics'] for i in path])
        return dfsGeneName

    def full_dfs_old(relation):
        source = [i[0] for i in relation]
        target = [i[1] for i in relation]
        rootIds = [i for i in source if i not in target]
        rootIds = list(set(rootIds))
        dfsPaths = []
        for root in rootIds:
            tempPaths = []
            relation_traverseFlag = {}
            for nodeArr in relation: relation_traverseFlag[Helper.edgeStr(nodeArr)] = False
            curPtr = root
            curPath = [curPtr]
            found, notEnd = False, False
            while len(curPath) != 0:
                if found: curPath.append(curPtr)
                found = False
                for nodeArr in relation:
                    if nodeArr[0] == curPtr:
                        if not relation_traverseFlag[Helper.edgeStr(nodeArr)] or Helper.pathNotTraveled(tempPaths, curPath):
                            curPtr = nodeArr[1]
                            relation_traverseFlag[Helper.edgeStr(nodeArr)] = True
                            found, notEnd = True, True
                            break
                if found: continue
                if notEnd: tempPaths.append(curPath.copy())
                notEnd = False
                curPath.pop()
                if len(curPath) != 0: curPtr = curPath[-1]
            dfsPaths += tempPaths
        return dfsPaths

    def source_dfs_old(relation, source):
        dfsPaths = []
        relation_traverseFlag = {}
        for nodeArr in relation: relation_traverseFlag[Helper.edgeStr(nodeArr)] = False
        curPtr = source
        curPath = [curPtr]
        found, notEnd = False, False
        while len(curPath) != 0:
            if found: curPath.append(curPtr)
            found = False
            for nodeArr in relation:
                if nodeArr[0] == curPtr:
                    if (not relation_traverseFlag[Helper.edgeStr(nodeArr)] or Helper.pathNotTraveled(dfsPaths, curPath)) and Helper.isNotCycle(curPath, nodeArr[1]):
                        curPtr = nodeArr[1]
                        relation_traverseFlag[Helper.edgeStr(nodeArr)] = True
                        found, notEnd = True, True
                        break
            if found: continue
            if notEnd: dfsPaths.append(curPath.copy())
            notEnd = False
            curPath.pop()
            if len(curPath) != 0: curPtr = curPath[-1]
        return dfsPaths

    def source_target_dfs_old(relation, source, target):
        sourceDFS = Traversal.source_dfs(relation, source)
        return [i for i in sourceDFS if i[-1] == target]

class PathCluster:
    def classify_pathways(dfsPaths):
        # classify pathways by disjoint set by their index
        entryClassifer = {}
        for i in range(0, len(dfsPaths)):
            for entry in dfsPaths[i]:
                if entry not in entryClassifer: entryClassifer[entry] = [i]
                else: entryClassifer[entry] += [i]
        entryClassifierList = [val for _, val in entryClassifer.items()]
        # merge the dfsPaths index by the node involved into the indexed pathway
        # growth algorithm: a list keep merging with other list, keep looping all lists until no further growth
        classifier = []
        for path in entryClassifierList:
            inclusive = False
            for C in classifier:
                if Helper.item_in_set(path, C):
                    inclusive = True
                    break
            if inclusive: continue

            start, growing = True, True
            while start or growing:
                growing, start = False, False
                for otherPath in entryClassifierList:
                    if path == otherPath: continue
                    if Helper.is_subset(otherPath, path): continue
                    if Helper.item_in_set(otherPath, path):
                        path = list(np.unique(np.asarray(path+otherPath)))
                        growing = True
            classifier.append(path)

        return classifier

    def classed_to_classic_path(classedPathsi, dfsPaths):
        classicPaths = []
        for C in classedPathsi:
            classicPaths.append([dfsPaths[entryID] for entryID in C])
        return classicPaths

    def to_relation(dfsPaths):
        relation = []
        for path in dfsPaths:
            for i in range(0, len(path)-1):
                if [path[i], path[i+1]] in relation: continue
                relation.append([path[i], path[i+1]])
        return relation
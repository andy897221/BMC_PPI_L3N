import bioGRID as bg
import os, json
import numpy as np

def parse_Giaever(root='./'):
    essentialGenes = []
    if os.path.exists(root+'data/parsed/giaever_essentialGenes.json'):
        with open(root+'data/parsed/giaever_essentialGenes.json', 'r') as f:
            essentialGenes = json.loads(f.read())
        return essentialGenes
    with open(root+'data/essential_gene/Giaever_Essential_ORFs.txt', 'r') as f:
        fLines = f.readlines()
        _, reverseGeneMap = bg.parse_geneName_map()
        for i in range(2, len(fLines)-4):
            sysGene = fLines[i].split('\t')[1].split(' ')[0]
            sysGeneVar = [sysGene.upper(), sysGene]
            for var in sysGeneVar:
                if var not in reverseGeneMap: continue
                essentialGenes.append(reverseGeneMap[var])
                break
    with open(root+'data/parsed/giaever_essentialGenes.json', 'w') as f:
            f.write(json.dumps(essentialGenes))
    return essentialGenes

def parse_YeastNet(root='./'):
    # is essentialPPI's both nodes are essentialGenes? or just one of them? need to confirm
    essentialPPI, essentialGenes = [], []
    if os.path.exists(root+'data/parsed/yeastNet_essentialGene.json'):
        with open(root+'data/parsed/yeastNet_essentialGenes.json', 'r') as f:
            essentialGenes = json.loads(f.read())
        return essentialGenes
    with open(root+'data/essential_gene/yeastNet_essentialGenes.txt', 'r') as f:
        _, reverseGeneMap = bg.parse_geneName_map()
        for line in f.readlines():
            line = line.split('\t')
            essentialPPI.append([reverseGeneMap[line[0]], reverseGeneMap[line[1]]])
    essentialGenes = list(set(np.asarray(essentialPPI).flatten()))
    with open(root+'data/parsed/yeastNet_essentialGenes.json', 'w') as f:
            f.write(json.dumps(essentialGenes))
    return essentialGenes

if __name__ == "__main__":
    essentialGenes = parse_Giaever()
    print(essentialGenes[-5:-1])
    essentialGenes = parse_YeastNet()
    print(essentialGenes[0:5])
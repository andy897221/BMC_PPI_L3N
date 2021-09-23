import json, os
import traversalHelper as tr

def parse_Pfam(filename="./data/Pfam/uniProt_Pfam_yeast_proteins.tab", wFile="./data/Pfam/parsed_geneName_Pfam.json", root="./"):
    # Pfam is domain ID
    # init index = 0, last col = Pfam, col 4 = gene name, skip first line
    # format: {geneName: [PfamID1, PfamID2, ...]}
    filename = root+filename
    wFile = root+wFile

    parsedPfam = {}
    if os.path.exists(wFile):
        with open(wFile) as f: parsedPfam = json.loads(f.read())
        return parsedPfam

    with open(filename, "r") as f:
        for line in f.readlines()[1:]:
            item = line[:-1].split("\t")
            geneName, PfamIDs = item[4].split(" ")[0], item[-1].split(";")[:-1]
            parsedPfam[geneName] = PfamIDs

    with open(wFile, "w") as f:
        f.write(json.dumps(parsedPfam))
    return parsedPfam

def parse_Pfam_DDI(ddiFile="./data/domine/INTERACTION.txt", wFile="./data/domine/parsed_domain.json", root="./"):
    # in ddiFile, HC can either be 'known interactions (PDB)' from db iPfam or 3did, or 'high confidence predictions'
    # MC means 'medium confidence", and LC means "low confidence"
    # ignore self-interaction
    # format: [[PfamID1, PfamID2], ...]
    ddiFile = root+ddiFile
    wFile = root+wFile

    DDIbr = []
    if os.path.exists(wFile):
        with open(wFile) as f: DDIbr = json.loads(f.read())
        return DDIbr

    with open(ddiFile, "r") as f:
        for line in f.readlines():
            line = line[:-1].split("|")
            if line[-2] != "HC": continue
            if line[0] == line[1]: continue
            DDIbr.append([line[0], line[1]])

    with open(wFile, "w") as f:
        f.write(json.dumps(DDIbr))
    return DDIbr

def _check_binaryDDI(p1, p2, GtoD, DDIr):
    if p1 not in GtoD or p2 not in GtoD: return False
    p1Domain, p2Domain = GtoD[p1], GtoD[p2]
    for p1d in p1Domain:
        if p1d not in DDIr: continue
        for p2d in p2Domain:
            if p2d in DDIr[p1d]: return True
    return False

def parse_binaryDDI_in_PPI(PPIs, GtoD_f='./data/Pfam/parsed_geneName_Pfam.json'
    , ddi_f='./data/domine/parsed_domain.json', root='./'):
    # determine if ddi exists for each ppi
    GtoD = parse_Pfam(wFile=GtoD_f, root=root)
    ddi = parse_Pfam_DDI(wFile=ddi_f, root=root)
    DDIr = tr.Helper.binary_to_relation(ddi, rSet=True)
    PPIr = tr.Helper.binary_to_relation(PPIs)
    ppi_DDIbinary = {}
    for protein in PPIr:
        ppi_DDIbinary[protein] = []
        for neighbor in PPIr[protein]:
            ppi_DDIbinary[protein].append(_check_binaryDDI(protein, neighbor, GtoD, DDIr))
    return PPIr, ppi_DDIbinary

if __name__ == "__main__":
    pass
    # parse_Pfam()
    # DDIbr = parse_Pfam_DDI()
    # print(DDIbr[0:10])
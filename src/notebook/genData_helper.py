import os
import json

def write_runTime(tag, runTime, fName="./linkPred_out/runTime.json"):
    if not os.path.exists(fName):
        with open(fName, 'w') as f: f.write(json.dumps({}))
    
    runTimeD = {}
    with open(fName, 'r') as f: runTimeD = json.loads(f.read())
    runTimeD[tag] = runTime
    with open(fName, 'w') as f: f.write(json.dumps(runTimeD))

def write_resultData(predPPI, predScore, fName, fLoc="./linkPred_out"):
    with open('{}/{}_PPI.json'.format(fLoc, fName), 'w') as f:
        f.write(json.dumps(predPPI))
    with open('{}/{}_score.json'.format(fLoc, fName), 'w') as f:
        f.write(json.dumps(predScore))
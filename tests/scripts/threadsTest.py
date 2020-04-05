"""Testing ccmap C extension

Usage:
    threadTest.py  cmap <threadNum> <dataSize> --lib <pathToPyLib> --pdbI <filepathFirstPdbFile> [ --pdbJ <filepathSecondPdbFile> ] [--encode --atomic --dist <ccDist>] 
    threadTest.py lcmap <threadNum> <dataSize> --lib <pathToPyLib> --pdbI <filepathFirstPdbFile> [ --pdbJ <filepathSecondPdbFile> ] [--encode --atomic --dist <ccDist>] 
    threadTest.py  zmap <threadNum> <dataSize> --lib <pathToPyLib> --inputs <transformationFile> --pdbI <filepathFirstPdbFile> --pdbJ <filepathSecondPdbFile> [--encode --atomic --dist <ccDist>]
    threadTest.py lzmap <threadNum> <dataSize> --lib <pathToPyLib> --inputs <transformationFile> --pdbI <filepathFirstPdbFile> --pdbJ <filepathSecondPdbFile> [--encode --atomic --dist <ccDist>]
    
    Options:
    -h --help
    <threadNum>  Number of threads (integer)
    <dataSize>  Number of structure to process (integer)
    --pdbI <filepathFirstPdbFile> Path to "receptor" PDB.
    --pdbJ <filepathSecondPdbFile> Path to "ligand" PDB.
    --inputs <transformationFile> Path to JSON transformation specs
    --lib <pathToPyLib>, folder containing the ccmap.so
    --encode  use integer encoding, default = False. Optional
    --atomic  compute atomic contact map, default = False. Optional
    --dist <ccDist> atomic contact threshold distance, in A.
"""

import sys, threading, json, time
from docopt import docopt

ARGS = docopt(__doc__, version="1.0.0")

sys.path.append(ARGS["--lib"])
import ccmap
import pyproteinsExt.structure.coordinates as PDB

print(ARGS)


parser = PDB.Parser()
pdbDictREC = parser.load(file=ARGS["--pdbI"]).atomDictorize
pdbDictLIG = parser.load(file=ARGS["--pdbJ"]).atomDictorize if ARGS["--pdbJ"] else None


vectors = None
if ARGS["--inputs"]:
    with open(ARGS["--inputs"], 'rb') as fp:
        vectors = json.load(fp)
        eulers     = [tuple(_) for _ in vectors['euler']]
        translations  = [tuple(_) for _ in vectors['translation']]
        ligOffset = vectors['ligOffset']
        recOffset = vectors['recOffset']
        print("Loaded Transformation to apply to \"", vectors['ligandFile'], "\"")

def cThread(*args, **kwargs):
    print(f"Starting  cThread {i}")
    print(f"{args} // f{kwargs}")
    tStart = time.time()
    #results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
    #        offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

def lcThread(*args, **kwargs):
    print(f"Starting lcThread {i}")
    print(f"{args} // f{kwargs}")
    tStart = time.time()
    #results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
    #        offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

def zThread(*args, **kwargs):
    print(f"Starting zThread {i}")
    print(f"{args} // f{kwargs}")
    tStart = time.time()
    #results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
    #        offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

def lzThread(*args, **kwargs):
    print(f"Starting lzThread {i}")
    print(f"{args} // f{kwargs}")
    #tStart = time.time()
    #results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
    #        offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

"""Yield successive interval boundaries spanning [0, iLen]"""
def splitInterval(iLen, nElem):
    assert(nElem <= iLen)
    nWidth = int(iLen/nElem)
    
    for i in range(nElem):
        top = (i+1) * nWidth
        if i == (nElem - 1): #and iLen%nElem != 0:
            top +=  iLen%nElem
        yield (i * nWidth, top)

def splitList(myList, nChunck):
    for x,y in splitInterval(len(myList), nChunck):
        yield myList[x:y]
    
threadNum = int(ARGS['<threadNum>'])
dataSize  = int(ARGS['<dataSize>'])

# Set worker function
wThread = None
if ARGS['cmap']:
    wThread = cThread 
elif ARGS['lcmap']:
    wThread = lcThread
elif ARGS['zmap']:
    wThread = zThread
else:
    wThread = lzThread
    
# Setting default thread named arguments
threadKwargs = {\
    "encode" : ARGS["--encode"],\
    "atomic" : ARGS["--atomic"] \
    } 
if ARGS['--dist']:
    threadKwargs["d"] = float(ARGS['--dist'])

# Setting default thread positional arguments
threadArgs = None

## We fill list-structures args with repetitions of the same PDBs
if ARGS['cmap'] or ARGS['lcmap']:
    wThread = cThread if ARGS['cmap'] else lcThread
    
    threadArgs = [[], []] if pdbDictLIG else [ [] ]
    cThread = 0
    for x,y in splitInterval(dataSize, threadNum):
        threadArgs[0].append([])
        if  pdbDictLIG:
            threadArgs[1].append([])
        for i in range(y - x):
            threadArgs[0][cThread].append(pdbDictREC)
            if  pdbDictLIG:
                threadArgs[1][cThread].append(pdbDictREC)
        cThread += 1

elif ARGS['zmap'] or ARGS['lzmap']:
    wThread = zThread if ARGS['zmap'] else lzThread
    kwargs = { 
    "offsetRec" : vectors['recOffset'], 
    "offsetLig" : vectors['ligOffset']
    }
    threadArgs = []
    for x,y in splitInterval(dataSize, threadNum):
        threadArgs.appends([ pdbDictREC, pdbDictLIG, \
            vectors['euler'][x:y],\
            vectors['translation'][x:y]\
        ])
    

mStart  = time.time()
dValues = [ 5.0 for i in range(threadNum) ]
output  = [ None for i in range(threadNum) ]
threadPool = [ threading.Thread(args = tuple(threadArgs[i]), kwargs = threadKwargs \
              , target=wThread) for i in range(threadNum)]

for th in threadPool:
    th.start()

for th in threadPool:
    th.join()

print(f"{threadNum} lz threads finished in { time.time() - mStart }")        

#if not bEncode:
#    for i,d in enumerate(output):
#        output[i] = json.loads(d)

#with open("threadsTest.json", 'w') as fp:
#    json.dump({ "threadData" : output }, fp)

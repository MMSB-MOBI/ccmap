"""Testing ccmap C extension

Usage:
    threadTest.py  cmap <threadNum> <dataSize> --lib <pathToPyLib> --pdbI <filepathFirstPdbFile> [ --pdbJ <filepathSecondPdbFile> ] [--encode --atomic --dist <ccDist> --out <fOut>] 
    threadTest.py lcmap <threadNum> <dataSize> --lib <pathToPyLib> --pdbI <filepathFirstPdbFile> [ --pdbJ <filepathSecondPdbFile> ] [--encode --atomic --dist <ccDist> --out <fOut>] 
    threadTest.py  zmap <threadNum> <dataSize> --lib <pathToPyLib> --inputs <transformationFile> --pdbI <filepathFirstPdbFile> --pdbJ <filepathSecondPdbFile> [--encode --atomic --dist <ccDist> --out <fOut>]
    threadTest.py lzmap <threadNum> <dataSize> --lib <pathToPyLib> --inputs <transformationFile> --pdbI <filepathFirstPdbFile> --pdbJ <filepathSecondPdbFile> [--encode --atomic --dist <ccDist> --out <fOut>]
    
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
    --out <fOut> path to json output file, default = \"threadsTest.json\"
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
    tArgs, tNum, results = args
    dual = len(tArgs) == 2
    results[tNum] = []
    print(f"Starting cThread {tNum} dual:{dual}")
    print( f"Shapes of ligand/receptor arrays {len(tArgs[0])}/{len(tArgs[1])}" )
    
    tStart = time.time()
    if dual:
        for rec,lig in zip(tArgs[0], tArgs[1]):
            results[tNum].append( ccmap.cmap(rec, lig, **kwargs) )
    else :
        for mol in tArgs[0]:
            results[tNum].append( ccmap.cmap(mol, **kwargs) )
    
    print(f"End of cThread {tNum} in { time.time() - tStart }")          
    return

def lcThread(*args, **kwargs): 
    tStart = time.time()
    tArgs, tNum, results = args
    assert len(tArgs) == 2
    

    results[tNum] = []
    print(f"Starting lcThread {tNum}")
    print( f"Shapes of ligand/receptor arrays {len(tArgs[0])}/{len(tArgs[1])}" )
    
    tStart = time.time()
    results[tNum] = ccmap.lcmap(*tArgs, **kwargs)
    
    print(f"End of lcThread {tNum} in { time.time() - tStart }")          
    return

def zThread(*args, **kwargs):
    tStart = time.time()
    tArgs, tNum, results = args
    pdbRec, pdbLib, eulerList, translationList = tArgs

    print(f"Starting zThread {tNum}, for {len(eulerList)} calls")
    
    results[tNum] = []
    for e,t in zip(eulerList, translationList):
        results[tNum].append( ccmap.zmap(pdbRec, pdbLib, e, t, **kwargs) )
        
    print(f"End of zThread {tNum} in { time.time() - tStart }")          
    return

def lzThread(*args, **kwargs):
    tStart = time.time()
    tArgs, tNum, results = args
    print(f"Starting lzThread {tNum}")
    
    results[tNum] = ccmap.lzmap(*tArgs, **kwargs)
        
    print(f"End of lzThread {tNum} in { time.time() - tStart }")          
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
    
    threadArgs = []
    # WRONG SHAPES [    [[] or ], .. NTHREADS]
    for x,y in splitInterval(dataSize, threadNum):
        threadArgs.append(( [], [] )) if pdbDictLIG else  threadArgs.append(( [] ))
        for i in range(y - x):
            threadArgs[-1][0].append(pdbDictREC)
            if pdbDictLIG:
                threadArgs[-1][1].append(pdbDictLIG)
        
elif ARGS['zmap'] or ARGS['lzmap']:
    if dataSize > len(vectors['euler']):
        print(f"dataSize {dataSize} exceeds avaible transformations {len(vectors['euler'])}, resizing it to maximum")
        dataSize = len(vectors['euler'])

    wThread = zThread if ARGS['zmap'] else lzThread
    threadKwargs["offsetRec"] = vectors['recOffset']
    threadKwargs["offsetLig"] = vectors['ligOffset']
   
    threadArgs = []
    for x,y in splitInterval(dataSize, threadNum):
        print(x,y)
        print(vectors['translation'][x:y])
        threadArgs.append( ( pdbDictREC, pdbDictLIG, \
                             vectors['euler'][x:y],\
                             vectors['translation'][x:y]\
                            ) )
    

mStart  = time.time()
output  = [ None for i in range(threadNum) ]
threadPool = [  threading.Thread(args = tuple( [ threadArgs[i], i, output ] )\
              , kwargs = threadKwargs \
              , target=wThread) \
            for i in range(threadNum)]

for th in threadPool:
    th.start()

for th in threadPool:
    th.join()

print(f"{threadNum} lz threads finished in { time.time() - mStart }")        

if not ARGS['--encode']:
    for i,d in enumerate(output):
        output[i] = json.loads(d)

fOut = ARGS['--out'] if ARGS['--out'] else "threadsTest.json"
with open(fOut, 'w') as fp:
    json.dump({ "threadData" : output }, fp)
  
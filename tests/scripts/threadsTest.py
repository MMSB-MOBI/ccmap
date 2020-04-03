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

print(f"Running {nThreads} threads of size {n} each")

parser = PDB.Parser()
pdbDictREC = parser.load(file=f"{ARGS["--pdbI"]}").atomDictorize
pdbDictLIG = parser.load(file=f"{ARGS["--pdbJ"]}").atomDictorize if ARGS["--pdbJ"] else None


vectors = None
if ARGS["--inputs"]:
    with open(ARGS["--inputs"], 'rb') as fp:
        vectors = json.load(fp)
        eulers     = [tuple(_) for _ in vectors['euler']]
        translations  = [tuple(_) for _ in vectors['translation']]
        ligOffset = vectors['ligOffset']
        recOffset = vectors['recOffset']
        print("Loaded Transformation to apply to \"", vectors['ligandFile'], "\"")

def cThread(d, pdb_lig, pdb_rec, e, t, lo, ro, results, i):
    print(f"Starting lz thread {i}")
    tStart = time.time()
    results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
            offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

def lcThread(d, pdb_lig, pdb_rec, e, t, lo, ro, results, i):
    print(f"Starting lz thread {i}")
    tStart = time.time()
    results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
            offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

def zThread(d, pdb_lig, pdb_rec, e, t, lo, ro, results, i):
    print(f"Starting lz thread {i}")
    tStart = time.time()
    results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
            offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

def lzThread(d, pdb_lig, pdb_rec, e, t, lo, ro, results, i):
    print(f"Starting lz thread {i}")
    tStart = time.time()
    results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
            offsetRec=ro, offsetLig=lo, distance=d, encode=bEncode)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

"""Yield successive interval boundaries spanning [0, iLen]"""
def splitInterval(iLen, nElem):
    assert(nElem <= iLen)
    
    iLen -= 1
    nWidth = int(iLen/nElem)
    
    for i in range(nElem):
        yield (i * nWidth, (i+1) * nWidth)
    
    if iLen%nElem != 0:
        yield (iLen - iLen%nElem , iLen + 1)
    
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
threadKwargs = [ {\
    "encode" : ARGS["--encode"],\
    "atomic" : ARGS["--atomic"] \
    } for t in range(threadNum) ]
for t in threadKwargs:
    if ARGS['--dist']:
        t["d"] = float(ARGS['--dist'])

# Setting default thread positional arguments
threadArgs = []
for t in threadKwargs:

if ARGS['cmap']:
    wThread = cThread 
    args = [ pdbDictREC ]
    if pdbDictLIG:
        args.append(pdbDictLIG)
elif ARGS['lcmap']:
    wThread = lcThread
    if pdbDictLIG:
        args = [ [ pdbDictREC for i in range(dataSize) ], \
                 [ pdbDictLIG for i in range(dataSize) ]  \
                ]
    else :
        args = [ [ pdbDictREC  for i in range(dataSize) ] ]

elif ARGS['zmap']:
    wThread = zThread
    args = [ pdbDictREC, pdbDictLIG, vectors['euler'][0], vectors['translation'] ]
else:
    wThread = lzThread
    args = [ pdbDictREC, pdbDictLIG ]
if ARGS["-zmap"] || ARGS["-lzmap"]:
    kwargs = { 
    "offsetRec" : vectors['recOffset'], 
    "offsetLig" : vectors['ligOffset']
    }
mStart  = time.time()
dValues = [ 5.0 for i in range(nThreads) ]
output  = [ None for i in range(nThreads) ]
threadPool = [threading.Thread(args=(d, pdbDictREC, pdbDictLIG,eulers[:n],\
              translations[:n], recOffset, ligOffset, output,i, ), target=lzThread) for i,d in enumerate(dValues)]

for th in threadPool:
    th.start()

for th in threadPool:
    th.join()

print(f"{nThreads} lz threads finished in { time.time() - mStart }")        

if not bEncode:
    for i,d in enumerate(output):
        output[i] = json.loads(d)

with open("threadsTest.json", 'w') as fp:
    json.dump({ "threadData" : output }, fp)
    ##fp.write("\"data\" : " + str(output) + "}")

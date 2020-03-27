import sys, threading, json, time
sys.path.append("/Users/guillaumelaunay/work/tmp/ccmap2/build/lib.macosx-10.9-x86_64-3.8")

import ccmap
import pyproteinsExt.structure.coordinates as PDB

nThreads = int(sys.argv[1])
n        = int(sys.argv[2])
print(f"Running {nThreads} threads of size {n} each")



folderRoot="/Users/guillaumelaunay/work/tmp/ccmap2/data/lzmap"

parser = PDB.Parser()
pdbDictREC = parser.load(file=f"{folderRoot}/1A2K_r_u.pdb").atomDictorize
pdbDictLIG = parser.load(file=f"{folderRoot}/1A2K_l_u.pdb").atomDictorize

with open(f"{folderRoot}/euler_translate_1A2K_1000.json", 'rb') as fp:
    vectors = json.load(fp)
eulers     = [tuple(_) for _ in vectors['euler']]
translations  = [tuple(_) for _ in vectors['translation']]

ligOffset = [-67.006, 0.11, -77.27]
recOffset = [-27.553, -8.229, -80.604]

def lzThread(d, pdb_lig, pdb_rec, e, t, lo, ro, results, i):
    print(f"Starting lz thread {i}")
    tStart = time.time()
    results[i] = ccmap.lzmap(pdb_rec, pdb_lig, e, t, \
            offsetRec=ro, offsetLig=lo, distance=d, encode=True)
    print(f"End of lz thread {i} in { time.time() - tStart }")          
    return

mStart  = time.time()
dValues = [ 5.0 for i in range(nThreads) ]
output  = [ None for i in range(nThreads) ]
threadPool = [threading.Thread(args=(d, pdbDictREC, pdbDictLIG,eulers[:n], translations[:n],ligOffset, recOffset, output,i, ), target=lzThread) for i,d in enumerate(dValues)]

for th in threadPool:
    th.start()

for th in threadPool:
    th.join()

print(f"{nThreads} lz threads finished in { time.time() - mStart }")        

with open("threadsTest.json", 'w') as fp:
    json.dump(output, fp)

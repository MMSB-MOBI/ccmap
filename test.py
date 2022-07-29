import sys
sys.path.insert(0,'/Users/guillaumelaunay/work/projects/ccmap/fibo_dev/ccmap/build/lib.macosx-10.9-x86_64-3.9')
import ccmap
from pypstruct import parseFilePDB

folderTests = "/Users/guillaumelaunay/work/projects/ccmap/fibo_dev/ccmap/tests/structures"

#pdbREC = parseFilePDB(filename=f"{folderTests}/1A2K_l_u.pdb")
pdbREC = parseFilePDB(filename=f"{folderTests}/small_peptide.pdb")
pdbDictREC = pdbREC.atomDictorize
noH_dict = {
    "x" : [],
    "y" : [],
    "z" : [],
    "seqRes" : [],
                    "chainID" : [],
                    "resName" : [],
                    "name" : []
}

for x,y,z,seqRes,chainID,resName, name in zip(pdbDictREC["x"], pdbDictREC["y"], pdbDictREC["z"],\
                                      pdbDictREC["seqRes"], pdbDictREC["chainID"],\
                                      pdbDictREC["resName"], pdbDictREC["name"]
                                    ):
    if name.startswith("H"):
        continue
    noH_dict["x"].append(x)
    noH_dict["y"].append(y)
    noH_dict["z"].append(z)
    noH_dict["seqRes"].append(seqRes)
    noH_dict["chainID"].append(chainID)    
    noH_dict["resName"].append(resName)
    noH_dict["name"].append(name)

#print(noH_dict["name"])

ccmap.sasa(noH_dict)

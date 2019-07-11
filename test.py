import sys, os
#sys.path.append("/Users/guillaumelaunay/work/DVL/pyproteinsExt/src")
#sys.path.append("/Users/guillaumelaunay/work/DVL/pyproteins/src")
#import pyproteinsExt.structure.coordinates as PDB
#import pyproteinsExt.structure.operations as PDBop

#parser = PDB.Parser()
#pdbObj = parser.load(file="/Users/guillaumelaunay/work/mesh/pdbFiles/2FJU.pdb")
#pdbObjRec = pdbObj.chain("A")
#pdbObjLig = pdbObj.chain("X")
import json
import ccmap


def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )

def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )

def _byteify(data, ignore_dicts = False):
    # if this is a unicode string, return its string representation
    if isinstance(data, unicode):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [ _byteify(item, ignore_dicts=True) for item in data ]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.iteritems()
        }
    # if it's anything else, return it in its original form
    return data





#atomDatumRec = None
#atomDatumLig = None
atomDatumRec = json_load_byteified(open('/home/glaunay/ccmap/data/2FJU_A.json'))
atomDatumLig = json_load_byteified(open('/home/glaunay/ccmap/data/2FJU_X.json'))
atomDatumRec_2 = json_load_byteified(open('/home/glaunay/ccmap/data/2FJU_A.json'))
atomDatumLig_2 = json_load_byteified(open('/home/glaunay/ccmap/data/2FJU_X.json'))
atomDatumRec_3 = json_load_byteified(open('/home/glaunay/ccmap/data/2FJU_A.json'))
atomDatumLig_3 = json_load_byteified(open('/home/glaunay/ccmap/data/2FJU_X.json'))
#with open('/home/glaunay/ccmap/data/2FJU_A.json') as json_data:
#    atomDatumRec = json.load(json_data)
#with open('/home/glaunay/ccmap/data/2FJU_X.json') as json_data:
#    atomDatumLig = json.load(json_data)

#print atomDatumLig
#jsonResults = ccmap.duals([(atomDatumRec, atomDatumLig),(atomDatumRec_2, atomDatumLig_2), (atomDatumRec_3, atomDatumLig_3)], 4.5)
#print jsonResults
jsonResults = ccmap.duals([(atomDatumRec, atomDatumLig),(atomDatumRec, atomDatumLig),(atomDatumRec, atomDatumLig),(atomDatumRec, atomDatumLig)], 4.5)
print sys.getrefcount(jsonResults)
print sys.getrefcount(atomDatumLig)
print sys.getrefcount(atomDatumRec)

#include "ccmapmodule_utils.h"
#include "mesh.h"
#include "atom_mapper.h"
#include "python_utils.h"
// --------------------- Utility Functions ---------------------  XREF SANITY ?

// COUNT REF NOT CHECKED
PyObject *buildPyValueSasaFrame(sasaFrame_t *sasaFrame){
    int len = sasaFrame->nbRes;
    PyObject *resnameList   = PyList_New(len);
    PyObject *resIDList     = PyList_New(len);
    PyObject *chainIDList   = PyList_New(len);
    
    int i = 0;
    residue_t *currResidue = sasaFrame->residueList->root;
    while(currResidue != NULL) {
        PyList_SetItem(resnameList,\
            i, Py_BuildValue("s", currResidue->resName));
        PyList_SetItem(resIDList,\
            i, Py_BuildValue("s", currResidue->resID));
        if (currResidue->ext_chainID != NULL)
            PyList_SetItem(chainIDList,\
                i, Py_BuildValue("s", currResidue->ext_chainID));
        else
            PyList_SetItem(chainIDList,\
                i, Py_BuildValue("C", (int)currResidue->chainID));
        currResidue = currResidue->nextResidueList;
        i++;
    }

    PyObject *sasaTupListFrame = PyList_New((Py_ssize_t) sasaFrame->nbFrame);
    PyObject *currSasaTupList  = NULL;
    PyObject *currSasaTup      = NULL;
    
    for (int i = 0 ; i < sasaFrame->nbFrame ; i++ ) {
        currSasaTupList = PyList_New(len); // DONT DECREF, i guess. -> Check its refcnt is one at the end of inner loop
        for (int j = 0 ; j < sasaFrame->nbRes ; j++ ) {
            currSasaTup = Py_BuildValue("(f,f)",\
                sasaFrame->sasa2upleArray[i][j].SASA,\
                sasaFrame->sasa2upleArray[i][j].frac);
            PyList_SetItem(currSasaTupList, j, currSasaTup); 
        }
        PyList_SetItem(sasaTupListFrame, i, currSasaTupList);
    }
    PyObject *results = PyDict_New();
    PyDict_SetItemString(results, "resname", resnameList); // Apparently doesnt steal ref of inserted obj, so we dont decref
    PyDict_SetItemString(results, "resID", resIDList);
    PyDict_SetItemString(results, "chainID", chainIDList);
    PyDict_SetItemString(results, "sasa", sasaTupListFrame);

    return results;
}

// Returns a representaion of ccmapVie as a PyList if bEncode, PyString otherwise
PyObject *ccmapViewToPyObject(ccmapView_t *ccmapView, bool bEncode) {
    
    if(bEncode) {
        PyObject *rValue = PyList_New(ccmapView->encodeLen);
        PyObject *pyValue;
        for (size_t i = 0 ; i < ccmapView->encodeLen ; i++){            
            pyValue = Py_BuildValue("i", ccmapView->asENCODE[i]);
            PyList_SetItem(rValue, i, pyValue );
        }
        return rValue;
    }
    return Py_BuildValue("s", ccmapView->asJSON);
}
// TO DO
/* pure C storage of current SASA calculations into sasaFrame accumulator */
void cmapViewAppendToSasaFrame(ccmapView_t *ccmapView, sasaFrame_t *sasaFrame, int iFrame){
    sasaResults_t *sasaResults = ccmapView->sasaResults;
    residue_sasa_t curResSasa;
    for (int i = 0; i < sasaResults->length ; i++) {
        curResSasa = sasaResults->residueSasaList[i];
        sasaFrame->sasa2upleArray[iFrame][i].SASA = curResSasa.nominal - curResSasa.buried;
        sasaFrame->sasa2upleArray[iFrame][i].frac = curResSasa.frac;
    }
}
//https://docs.python.org/3/c-api/arg.html#building-values
// https://web.mit.edu/people/amliu/vrut/python/ext/buildValue.html
PyObject *ccmapViewToSasaPyDict(ccmapView_t *ccmapView) {
    sasaResults_t *sasaResults = ccmapView->sasaResults;
    residue_sasa_t curResSasa;
    PyObject *pyValue;
    char chainID[2];
    chainID[1] = '\0';
    Py_ssize_t len = sasaResults->length;
    PyObject *mainDict = PyDict_New();
    const char* mKey = "freeASA"; 
    PyObject *mainList = PyList_New(len);
    PyDict_SetItemString( mainDict, mKey , mainList);
    for (int i = 0; i < sasaResults->length ; i++) {
        curResSasa = sasaResults->residueSasaList[i];
        chainID[0] = curResSasa.chainID;

        pyValue = Py_BuildValue("{s:s,s:s,s:s,s:f,s:f}",
                      "resname" , curResSasa.resname,
                       "resID"  , curResSasa.resID, 
                       "chainID", chainID,
                       "SASA"   , curResSasa.nominal - curResSasa.buried,
                       "frac"   ,  curResSasa.frac );
        PyList_SetItem(mainList, i, pyValue); 
    }
//    fprintf(stderr, "mainDict:%zd\n", Py_REFCNT(mainDict) );
 //   fprintf(stderr, "mainList:%zd\n", Py_REFCNT(mainList) ); //- > 2
    Py_DECREF(mainList);
   
    return mainDict;
}

// NEED TO CHECK XREFs !!!
PyObject *ccmapViewsToPyObject(ccmapView_t **ccmapViews, int nViews, bool bEncode) {
    PyObject *rValue;
    ccmapView_t *currView;
    if(bEncode) {
        #ifdef DEBUG
        PySys_WriteStdout("Building encoding results\n");
        #endif
        rValue = PyList_New(nViews);
        PyObject *pyCurrList, *pyValue;
        
        for (int iView = 0 ; iView < nViews ; iView++) {
            currView = ccmapViews[iView];
            PyList_SetItem( rValue, iView, PyList_New(currView->encodeLen) );
            pyCurrList = PyList_GetItem(rValue, iView); // Stolen
            for (size_t i = 0 ; i < currView->encodeLen ; i++){            
                pyValue = Py_BuildValue("i", currView->asENCODE[i]);
                PyList_SetItem(pyCurrList, i, pyValue); //pyValue Ref(v=1) is stolen by pyCurrList
            }
        }
        return rValue;
    }
   
    string_t *jsonStringEncodeMany = createString();
    jsonStringEncodeMany->append(jsonStringEncodeMany, "{\"type\":\"lzmap\", \"data\":[");
    for (int iView = 0 ; iView < nViews ; iView++) {
            jsonStringEncodeMany->append(jsonStringEncodeMany, ccmapViews[iView]->asJSON);
        if (iView < nViews - 1)
            jsonStringEncodeMany->append(jsonStringEncodeMany, ",");  
    }
    jsonStringEncodeMany->append(jsonStringEncodeMany, "]}");

    rValue = Py_BuildValue("s", jsonStringEncodeMany->value);
    destroyString(jsonStringEncodeMany);
    return rValue;
}

int backMapCoordinates(atom_t *atomListRoot,  PyObject *pyDictObject) {
    PyObject *pItem;
    atom_t *atomCurrent = atomListRoot;
    PyObject* pyObj_x = PyDict_GetItemString(pyDictObject, "x");
    PyObject* pyObj_y = PyDict_GetItemString(pyDictObject, "y");
    PyObject* pyObj_z = PyDict_GetItemString(pyDictObject, "z");
    Py_ssize_t n = PyList_Size(pyObj_x);

    int atomIndex = 0;

    double test=-999.9;

     while (atomCurrent != NULL) {
        pItem = Py_BuildValue("d", atomCurrent->x);// New Ref
        PyObject_AsDouble(pItem, &test);
        PyList_SetItem(pyObj_x, atomIndex, pItem); // Stolen

        pItem = Py_BuildValue("d", atomCurrent->y);
        PyObject_AsDouble(pItem, &test);
        PyList_SetItem(pyObj_y, atomIndex, pItem);

        pItem = Py_BuildValue("d", atomCurrent->z);
        PyObject_AsDouble(pItem, &test);
        PyList_SetItem(pyObj_z, atomIndex, pItem);

        atomCurrent = atomCurrent->nextAtomList;
        atomIndex++;
    }

    return 1;
}

void freeBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n) {
    PyMem_Free(x);
    PyMem_Free(y);
    PyMem_Free(z);
    PyMem_Free(chainID);
    for (int i = 0; i < n ; i++) {
        PyMem_Free(resID[i]);
        PyMem_Free(resName[i]);
        PyMem_Free(name[i]);
    }
    PyMem_Free(resID);
    PyMem_Free(resName);
    PyMem_Free(name);
}

// Load content of a Python dictionary of residue atom radii into a atomMapper
// Expected to be loaded once per "thread"
atom_map_t *dictRadiiToAtomMapper(PyObject *pyRadiiDictObject) {
    PyObject *key, *value;
    PyObject *pItem_atom_name_radius_tuple;
    Py_ssize_t pos = 0;

    char resname[81];
    int i_atom, nb_atom;
    PyObject *pyObj_atom_name, *pyObj_atom_radius;
    
    char **atom_name_buffer;
    float *atom_radii_buffer;
    create_buffers(&atom_name_buffer, &atom_radii_buffer);
    atom_map_t *atom_map = createAtomMapper();
   // fprintf(stderr, "Original pyRadiiDictObject:%zd\n", Py_REFCNT(pyRadiiDictObject) );
            
    // Atom data as list or tuple
    PyObject *(*get_item) (PyObject *, Py_ssize_t) = NULL;
    while (PyDict_Next(pyRadiiDictObject, &pos, &key, &value)) {
        //fprintf(stderr, "key:%zd\n", Py_REFCNT(key) );
        Py_INCREF(key);
        PyObject_ToChar(key, resname);
        Py_DECREF(key);
        //fprintf(stderr, "key:%zd\n", Py_REFCNT(key) );
        nb_atom = PyList_Size(value);     
        #ifdef DEBUG
        PySys_WriteStderr("Current Key is %s [%d atoms]\n", resname, nb_atom);
        #endif
        for (i_atom = 0 ; i_atom < nb_atom ; i_atom++) {
            pItem_atom_name_radius_tuple = PyList_GetItem(value, i_atom);
           // fprintf(stderr, "atom_name_radius_tuple:%zd\n", Py_REFCNT(pItem_atom_name_radius_tuple) );
            Py_INCREF(pItem_atom_name_radius_tuple);
          //  fprintf(stderr, "atom_name_radius_tuple:%zd\n", Py_REFCNT(pItem_atom_name_radius_tuple) );
            if(get_item == NULL)
                get_item = PyTuple_Check(pItem_atom_name_radius_tuple) ? \
                            PyTuple_GetItem :\
                            PyList_GetItem;


            pyObj_atom_name   = get_item(pItem_atom_name_radius_tuple, 0); // Py_INCREF // Py_DECREF ??
            Py_INCREF(pyObj_atom_name);
            PyObject_ToChar(pyObj_atom_name, atom_name_buffer[i_atom]);
            pyObj_atom_radius = get_item(pItem_atom_name_radius_tuple, 1);
            Py_INCREF(pyObj_atom_radius);
            atom_radii_buffer[i_atom] =  (float) PyFloat_AsDouble(pyObj_atom_radius);
            if(atom_radii_buffer[i_atom] > atom_map->maxRadius)
                atom_map->maxRadius = atom_radii_buffer[i_atom];
            Py_DECREF(pyObj_atom_name);
            Py_DECREF(pyObj_atom_radius);       

            Py_DECREF(pItem_atom_name_radius_tuple);
            
          /*
           fprintf(stderr,\
            "atom loop members:\n\tpyObj_atom_name:%zd\n\pyObj_atom_radius:%zd\n\pItem_atom_name_radius_tuple:%zd\n",\
            Py_REFCNT(pyObj_atom_name),  Py_REFCNT(pyObj_atom_radius),  Py_REFCNT(pyObj_atom_radius) );
            */
        }
        addMapGroup(atom_map, atom_name_buffer, resname, atom_radii_buffer, nb_atom);
       
    }
    //Py_DECREF(pyRadiiDictObject); // Do not do that !!
    //fprintf(stderr, "==>pyRadiiDictObject:%zd\n", Py_REFCNT(pyRadiiDictObject) );
    destroy_buffers(atom_name_buffer, atom_radii_buffer);
    return atom_map;
}

// Convert a Python dictionnarized structure into a list of atoms, expecting no hydrogens.
// bASA Shall be a optional atomMapper pointer
atom_t *structDictToAtoms(PyObject *pyDictObject, int *nAtoms, float probeRadius, atom_map_t *aMap) {
#ifdef DEBUG
    char DBG_buffer[1024];
    sprintf(DBG_buffer, "Running structDictToAtoms: bASA:%s, probeRadius:%f\n", aMap != NULL ? "true": "false", probeRadius);
     
    #ifdef AS_PYTHON_EXTENSION
        PySys_WriteStderr("%s", DBG_buffer);
    #elif
        fprintf(stderr, "%s", DBG_buffer);
    #endif
#endif


#ifdef PYMEM_CHECK
    fprintf(stderr, "structDictToAtoms entry MEMORY SUMMARY\n");
    fprintf(stderr, "pyDictObject:%zd\n", Py_REFCNT(pyDictObject) );
#endif
    
    PyObject* pyObj_x        = PyDict_GetItemString(pyDictObject, "x");
    Py_INCREF(pyObj_x);
    PyObject* pyObj_y        = PyDict_GetItemString(pyDictObject, "y");
    Py_INCREF(pyObj_y);
    PyObject* pyObj_z        = PyDict_GetItemString(pyDictObject, "z");
    Py_INCREF(pyObj_z);
    PyObject* pyObj_chainID  = PyDict_GetItemString(pyDictObject, "chainID");
    Py_INCREF(pyObj_chainID);
    PyObject* pyObj_resSeq   = PyDict_GetItemString(pyDictObject, "seqRes");
    Py_INCREF(pyObj_resSeq);
    PyObject* pyObj_resName  = PyDict_GetItemString(pyDictObject, "resName");
    Py_INCREF(pyObj_resName);
    PyObject* pyObj_atomName = PyDict_GetItemString(pyDictObject, "name");
    Py_INCREF(pyObj_atomName);
    Py_ssize_t n             = PyList_Size(pyObj_x);
    *nAtoms = (int) n;

    /*
    All unpackXX calls do memory allocation, which needs subsequent common call to freeBuffer()
    */
 
    double *coorX = unpackCoordinates(pyObj_x);
    Py_DECREF(                        pyObj_x);
    double *coorY = unpackCoordinates(pyObj_y);
    Py_DECREF(                        pyObj_y);
    double *coorZ = unpackCoordinates(pyObj_z);
    Py_DECREF(                       pyObj_z); 
    
    char *chainID;
    unpackChainID(pyObj_chainID, &chainID);
    Py_DECREF(    pyObj_chainID);
   
    char **resSeq;
    unpackString(pyObj_resSeq, &resSeq);
    Py_DECREF(   pyObj_resSeq);
   
    char **resName;
    unpackString(pyObj_resName, &resName);
    Py_DECREF(   pyObj_resName);
   
    char **atomName;
    unpackString(pyObj_atomName, &atomName);
    Py_DECREF(   pyObj_atomName);

    #ifdef DEBUG
        sprintf(DBG_buffer, "structDictToAtoms: Calling readFromArrays over %d atoms\n", *nAtoms);
        printOnContextStderr(DBG_buffer);
    #endif
    // Safe here
    /* Create data structures and compute */
    atom_t *atomList = readFromArrays(*nAtoms, coorX, coorY, coorZ, chainID, resSeq, resName, atomName, aMap, probeRadius);
    
    freeBuffers(coorX, coorY, coorZ, chainID, resSeq, resName,  atomName, *nAtoms);

    
#ifdef PYMEM_CHECK
    fprintf(stderr, "structDictToAtoms exit MEMORY SUMMARY\n");
    fprintf(stderr, "pyObj_x:%zd\n", Py_REFCNT(pyObj_x) );
    fprintf(stderr, "pyObj_y:%zd\n", Py_REFCNT(pyObj_y) );
    fprintf(stderr, "pyObj_z:%zd\n", Py_REFCNT(pyObj_z) );
    fprintf(stderr, "pyObj_chainID:%zd\n", Py_REFCNT(pyObj_chainID) );
    fprintf(stderr, "pyObj_resSeq:%zd\n", Py_REFCNT(pyObj_resSeq) );
    fprintf(stderr, "pyObj_resName:%zd\n", Py_REFCNT(pyObj_resName) );
    fprintf(stderr, "pyObj_atomName:%zd\n", Py_REFCNT(pyObj_atomName) );
    fprintf(stderr, "pyDictObject:%zd\n", Py_REFCNT(pyDictObject));
#endif

    #ifdef DEBUG
        printOnContextStderr("structDictToAtoms: All buffer free, exiting\n");
    #endif

    return atomList;
}



void setBooleanFromParsing(PyObject *pyObjectBool, bool *bResults) 
{
*bResults = false;
if (pyObjectBool != NULL) {
    if (PyObject_IsTrue(pyObjectBool))
        *bResults = true;
    //Py_DECREF(pyObjectBool); It's borrowed from main fun caller and we dont need it anymore, DONT TOUCH IT!
    }
#ifdef PYMEM_CHECK
    if (pyObjectBool != NULL)
        fprintf(stderr, "setBooleanFromParsing pyObjectBool refcount %zd\n", pyObjectBool->ob_refcnt);
#endif
}


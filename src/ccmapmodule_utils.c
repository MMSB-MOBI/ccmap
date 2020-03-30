#include "ccmapmodule_utils.h"
#include "mesh.h"

// --------------------- Utility Functions ---------------------  XREF SANITY ?

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

// NEED TO CHECK XREFs !!!
PyObject *ccmapViewsToPyObject(ccmapView_t **ccmapViews, int nViews, bool bEncode) {
    PyObject *rValue;
    ccmapView_t *currView;
    if(bEncode) {
        PySys_WriteStdout("Building encoding results\n");
        rValue = PyList_New(nViews);
        PyObject *pyCurrList, *pyValue;
        
        for (int iView = 0 ; iView < nViews ; iView++) {
            currView = ccmapViews[iView];
           // PySys_WriteStdout("CurrView is of length %lu\n", currView->encodeLen);
            PyList_SetItem( rValue, iView, PyList_New(currView->encodeLen) );
            pyCurrList = PyList_GetItem(rValue, iView);
            for (size_t i = 0 ; i < currView->encodeLen ; i++){            
                pyValue = Py_BuildValue("i", currView->asENCODE[i]); // pyValue xref ?
                PyList_SetItem(pyCurrList, i, pyValue);
            }
        }
        return rValue;
    }
    //return Py_BuildValue("s", "RTOTO");
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

int PyObject_AsDouble(PyObject *py_obj, double *x)
{
  
  PyObject *py_float = PyNumber_Float(py_obj);

  if (py_float == NULL) return -1;
  *x = 0; 
  *x += PyFloat_AsDouble(py_float);

 //*x = 5.0;
  Py_DECREF(py_float);
  //PySys_WriteStdout("REF COUNT results  :: is %d\n", Py_REFCNT(py_float) ); # IT IS STILL 1
  return 0;
}

int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size) {
    int i;

    if (py_list == NULL) return 1;

    if (!PyList_Check(py_list)) return 1;

    if (size != PyList_Size(py_list)) return 1;

    for (i=0; i<size; i++) {
        PyObject *py_float = PyList_GetItem(py_list, i);
        if (py_float == NULL || PyObject_AsDouble(py_float, &(x[i])))
        return 1;
    }

    return 0;
}

int backMapCoordinates(atom_t *atomListRoot,  PyObject *pyDictObject) {
    PyObject *pItem;
    Py_ssize_t n;
    atom_t *atomCurrent = atomListRoot;
    PyObject* pyObj_x = PyDict_GetItemString(pyDictObject, "x");
    PyObject* pyObj_y = PyDict_GetItemString(pyDictObject, "y");
    PyObject* pyObj_z = PyDict_GetItemString(pyDictObject, "z");
    n = PyList_Size(pyObj_x);

    int atomIndex = 0;

    double test=-999.9;

     while (atomCurrent != NULL) {
        pItem = Py_BuildValue("d", atomCurrent->x);
        PyObject_AsDouble(pItem, &test);
        PyList_SetItem(pyObj_x, atomIndex, pItem);

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

double **createListVector3(PyObject *pyObject_List, Py_ssize_t *len) {
    *len = PyList_Size(pyObject_List);  
    double **vList = PyMem_New(double*, *len);
    
    PyObject *currPyTuple;
    bool bError =false;
    for (int i = 0 ; i < *len ; i++) {
        currPyTuple = PyList_GetItem(pyObject_List, i);
        vList[i] = PyMem_New(double, 3);
        vList[i] = unpackVector3(currPyTuple);
        if (vList[i] == NULL) {
            bError = true;
            break;
        }
    }
    if (bError)
        vList = destroyListVector3(vList, *len);
    return vList;
}

double **destroyListVector3(double **vList, Py_ssize_t len) {
    for (int i = 0 ; i < len ; i++)
        destroyVector3(vList[i]);
    PyMem_Free(vList);
    return vList;
}

double *unpackVector3(PyObject *pyObject) {
#ifdef DEBUG
    PySys_WriteStdout("--->Unpack Vector3\n");
#endif
    PyObject *(*PyArray_GetItem)(PyObject *, Py_ssize_t);
    Py_ssize_t (*PyArray_Size)(PyObject *);

    if(PyList_Check(pyObject)) {
        PyArray_GetItem = &PyList_GetItem;
        PyArray_Size = &PyList_Size;
    } else if(PyTuple_Check(pyObject)) {
        PyArray_GetItem = &PyTuple_GetItem;
        PyArray_Size = &PyTuple_Size;
    } else {
        PyErr_SetString(PyExc_TypeError, "Error unpacking a vector3 from unknown array type");
        return NULL;
    }
    if (PyArray_Size(pyObject) != 3) {
        PyErr_SetString(PyExc_TypeError, "Error unpacking a vector3 from a python array of size != 3");
        return NULL;
    }

    double *vector = PyMem_New(double, 3);
    PyObject *pItem;
   
    for (int i = 0 ; i < 3 ; i++) {
        pItem = PyArray_GetItem(pyObject, i);
        if(!PyFloat_Check(pItem)) {
            PyMem_Free(vector);
            PyErr_SetString(PyExc_TypeError, "3D vector element items must be float.");
            return NULL;
        }
        PyObject_AsDouble( pItem, &(vector[i]));
    }
    return vector;
}

double *destroyVector3(double *vector) {
    if (vector == NULL)
        return vector;
    PyMem_Free(vector);
    return vector;
}

int unpackChainID(PyObject *pListChainID, char **buffer) {
#ifdef DEBUG
    PySys_WriteStdout("--->Unpack chainID\n");
#endif
    PyObject *pItem;
    Py_ssize_t n;
    int i;

    n = PyList_Size(pListChainID);
   
    *buffer = PyMem_New(char, n);

    PyObject* objectsRepresentation;
    const char* s;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListChainID, i);

        objectsRepresentation = PyObject_Repr(pItem); //Now a unicode object
        PyObject* pyStr = PyUnicode_AsUTF8String(objectsRepresentation);
        s = PyBytes_AS_STRING(pyStr);
        Py_XDECREF(pyStr);

        (*buffer)[i] = s[1];
    }
    Py_XDECREF(pItem); // LAST MOD
    Py_XDECREF(objectsRepresentation);
    return 1;
}

int unpackString(PyObject *pListOfStrings, char ***buffer) {
#ifdef DEBUG
    PySys_WriteStdout("--->Unpack string\n");
#endif
    PyObject *pItem;
    Py_ssize_t n;
    int i;

    //char **buffer = *_buffer;
    n = PyList_Size(pListOfStrings);
    *buffer = PyMem_New(char*, n);

    PyObject* objectsRepresentation;
    const char* s;
    int sLen;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListOfStrings, i);

        objectsRepresentation = PyObject_Repr(pItem); //Now a unicode object
        PyObject* pyStr = PyUnicode_AsUTF8String(objectsRepresentation);
        s = PyBytes_AS_STRING(pyStr);
        Py_XDECREF(pyStr);

        sLen =  strlen(s); // This corresponds to the actual string surrounded by \' , ie : 'MYTSRING'
        //PySys_WriteStdout("--->%s[%d]\n", s, strlen(s));
        (*buffer)[i] = PyMem_New(char, sLen - 1);
        for (int j = 1 ; j < sLen - 1 ; j++) {
            (*buffer)[i][j - 1] = s[j];
        }
        (*buffer)[i][sLen - 2] = '\0';
        Py_XDECREF(objectsRepresentation);
        //PySys_WriteStderr("NO DECFREF");
        //PySys_WriteStdout("translated to --->\"%s\"[%d]\n", (*buffer)[i], sLen - 1);
       // PySys_WriteStdout("translated to --->%s[%d]\n", (*buffer)[i]);
    }
#ifdef DEBUG
    PySys_WriteStdout("REF COUNT :: is %d\n", Py_REFCNT(objectsRepresentation) );
#endif
    return 1;
}

double *unpackCoordinates(PyObject *pListCoor) {
    PyObject *pItem;
    Py_ssize_t n = PyList_Size(pListCoor);
    double *buffer = PyMem_New(double, n);

    for (int i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListCoor, i);
        if(!PyFloat_Check(pItem)) {
            PyErr_SetString(PyExc_TypeError, "coordinate items must be float.");
            PyMem_Free(buffer);
            return NULL;
        }

        PyObject_AsDouble(pItem, &(buffer[i]) );
        //PySys_WriteStdout("TEST:: %.2f\n", buffer[i] );
    }
    #ifdef DEBUG
    PySys_WriteStderr("Allocation done\n");
    #endif
    return buffer;
}

void freeBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n) {
#ifdef DEBUG
    PySys_WriteStdout("Freeing all I buffers of size %d\n", n);
#endif
    //fprintf(stderr, "Freeing all I buffers of size %d\n", n);
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
#ifdef DEBUG
    PySys_WriteStdout("Done\n");
#endif
}

atom_t *structDictToAtoms(PyObject *pyDictObject, int *nAtoms) {
#ifdef DEBUG
    PySys_WriteStdout("\n\n=======================\n");
    PySys_WriteStdout("REF COUNT current Dict :: is %d\n", Py_REFCNT(pyDictObject));
#endif
    //Return value: Borrowed reference.
    PyObject* pyObj_x = PyDict_GetItemString(pyDictObject, "x");
    PyObject* pyObj_y = PyDict_GetItemString(pyDictObject, "y");
    PyObject* pyObj_z = PyDict_GetItemString(pyDictObject, "z");
    PyObject* pyObj_chainID = PyDict_GetItemString(pyDictObject, "chainID");
    PyObject* pyObj_resSeq = PyDict_GetItemString(pyDictObject, "seqRes");
    PyObject* pyObj_resName = PyDict_GetItemString(pyDictObject, "resName");
    PyObject* pyObj_atomName = PyDict_GetItemString(pyDictObject, "name");

    Py_ssize_t n = PyList_Size(pyObj_x);
    *nAtoms = (int) n;
#ifdef DEBUG
    PySys_WriteStdout("Unpacking a %d atoms structure dictionary\n", *nAtoms);
#endif
    /*
    All unpackXX calls do memory allocation, which needs subsequent common call to freeBuffer()
    */
    double *coorX = unpackCoordinates(pyObj_x);
    double *coorY = unpackCoordinates(pyObj_y);
    double *coorZ = unpackCoordinates(pyObj_z);
    /* DONT DECREF REFERENCE IS BORROWED !
    Py_DECREF(pyObj_x);
    Py_DECREF(pyObj_y);
    Py_DECREF(pyObj_z);
    */
    char *chainID;
    unpackChainID(pyObj_chainID, &chainID);
    //Py_DECREF(pyObj_chainID);

    char **resSeq;
    unpackString(pyObj_resSeq, &resSeq);
    //Py_DECREF(pyObj_resSeq);

    char **resName;
    unpackString(pyObj_resName, &resName);
    //Py_DECREF(pyObj_resName);

    char **atomName;
    unpackString(pyObj_atomName, &atomName);
    //Py_DECREF(pyObj_atomName);

    /* Create data structures and compute */
    atom_t *atomList = readFromArrays(*nAtoms, coorX, coorY, coorZ, chainID, resSeq, resName, atomName);

    freeBuffers(coorX, coorY, coorZ, chainID, resSeq, resName,  atomName, *nAtoms);

#ifdef DEBUG
    PySys_WriteStdout("REF COUNT X :: is %d\n", Py_REFCNT(pyObj_x) );
    PySys_WriteStdout("REF COUNT Y :: is %d\n", Py_REFCNT(pyObj_y) );
    PySys_WriteStdout("REF COUNT Z :: is %d\n", Py_REFCNT(pyObj_z) );
    PySys_WriteStdout("REF COUNT chainID :: is %d\n", Py_REFCNT(pyObj_chainID) );
    PySys_WriteStdout("REF COUNT resSeq :: is %d\n", Py_REFCNT(pyObj_resSeq) );
    PySys_WriteStdout("REF COUNT resName :: is %d\n", Py_REFCNT(pyObj_resName) );
    PySys_WriteStdout("REF COUNT atomName :: is %d\n", Py_REFCNT(pyObj_atomName) );

    PySys_WriteStdout("REF COUNT current Dict :: is %d\n", Py_REFCNT(pyDictObject));
    PySys_WriteStdout("Returning atomList\n");
#endif
    return atomList;
}



void setBooleanFromParsing(PyObject *pyObjectBool, bool *bResults) 
{
*bResults = false;
if (pyObjectBool != NULL) {
    if (PyObject_IsTrue(pyObjectBool))
        *bResults = true;
    Py_XDECREF(pyObjectBool);
    }
}

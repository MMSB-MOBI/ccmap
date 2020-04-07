#include "ccmapmodule_utils.h"
#include "mesh.h"

// --------------------- Utility Functions ---------------------  XREF SANITY ?

bool PyArray_Equal(PyObject *arrayI, PyObject *arrayJ) {
    Py_ssize_t (*PyArray_SizeI)(PyObject *);
    Py_ssize_t (*PyArray_SizeJ)(PyObject *);

    PyArray_SizeI = PyList_Check(arrayI) ? &PyList_Size : &PyTuple_Size;
    PyArray_SizeJ = PyList_Check(arrayJ) ? &PyList_Size : &PyTuple_Size;
  
    return PyArray_SizeI == PyArray_SizeJ;
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

int PyObject_AsDouble(PyObject *py_obj, double *x)
{
  
  PyObject *py_float = PyNumber_Float(py_obj); // New Ref

  if (py_float == NULL) return -1;
  *x = 0; 
  *x += PyFloat_AsDouble(py_float);

  Py_DECREF(py_float);
 
  return 0;
}

int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size) {
    int i;

    if (py_list == NULL) return 1;

    if (!PyList_Check(py_list)) return 1;

    if (size != PyList_Size(py_list)) return 1;

    for ( i = 0 ; i < size ; i++ ) {
        PyObject *py_float = PyList_GetItem(py_list, i); // Borrowed Ref
        Py_XINCREF(py_float);
        if (py_float == NULL || PyObject_AsDouble(py_float, &(x[i]))) 
            return 1;
        Py_DECREF(py_float);
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
/* Get item at position in a tuple or a list */
PyObject *PyArray_GetItem(PyObject *pyObject_array, Py_ssize_t position) {
    PyObject *(*_PyArray_GetItem)(PyObject *, Py_ssize_t);

    if(PyList_Check(pyObject_array)) {
        _PyArray_GetItem = &PyList_GetItem;
    } else if(PyTuple_Check(pyObject_array)) {
        _PyArray_GetItem = &PyTuple_GetItem;
    } else {
        return NULL;
    }
    return _PyArray_GetItem(pyObject_array, position);
}  

bool PyArray_Check(PyObject *pyObject_arrayMaybe) {
    if (PyList_Check(pyObject_arrayMaybe) )
        return true;
    if (PyTuple_Check(pyObject_arrayMaybe) )
        return true;
    
    return false;
}

Py_ssize_t PyArray_Size(PyObject *pyObject_array) {
    Py_ssize_t (*_PyArray_Size)(PyObject *);

    if(PyList_Check(pyObject_array)) {
        _PyArray_Size = &PyList_Size;
    } else if(PyTuple_Check(pyObject_array)) {
        _PyArray_Size = &PyTuple_Size;
    } else {
        return 0; // ??
    }
    return _PyArray_Size(pyObject_array);
}

double **createListVector3(PyObject *pyObject_array, Py_ssize_t *len) {
    PyObject *(*_PyArray_GetItem)(PyObject *, Py_ssize_t);
    Py_ssize_t (*_PyArray_Size)(PyObject *);

    if(PyList_Check(pyObject_array)) {
        _PyArray_GetItem = &PyList_GetItem;
        _PyArray_Size = &PyList_Size;
    } else if(PyTuple_Check(pyObject_array)) {
        _PyArray_GetItem = &PyTuple_GetItem;
        _PyArray_Size = &PyTuple_Size;
    } else {
        PyErr_SetString(PyExc_TypeError, "Error unpacking list of vector3 from unknown array type");
        return NULL;
    }
    *len = _PyArray_Size(pyObject_array);  
    double **vList = PyMem_New(double*, *len);
    
    PyObject *currPyArray;
    bool bError = false;
    for (int i = 0 ; i < *len ; i++) {
        currPyArray = _PyArray_GetItem(pyObject_array, i);
        Py_INCREF(currPyArray);
        
        vList[i] = PyMem_New(double, 3);
        vList[i] = unpackVector3(currPyArray);
        
        Py_DECREF(currPyArray);

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

    PyObject *(*_PyArray_GetItem)(PyObject *, Py_ssize_t);
    Py_ssize_t (*_PyArray_Size)(PyObject *);

    if(PyList_Check(pyObject)) {
        _PyArray_GetItem = &PyList_GetItem;
        _PyArray_Size = &PyList_Size;
    } else if(PyTuple_Check(pyObject)) {
        _PyArray_GetItem = &PyTuple_GetItem;
        _PyArray_Size = &PyTuple_Size;
    } else {
        PyErr_SetString(PyExc_TypeError, "Error unpacking a vector3 from unknown array type");
        return NULL;
    }
    if (_PyArray_Size(pyObject) != 3) {
        PyErr_SetString(PyExc_TypeError, "Error unpacking a vector3 from a python array of size != 3");
        return NULL;
    }

    double *vector = PyMem_New(double, 3);
    PyObject *pItem;
   
    for (int i = 0 ; i < 3 ; i++) {
        pItem = _PyArray_GetItem(pyObject, i);
        Py_INCREF(pItem);
        if(!PyFloat_Check(pItem)) {
            PyMem_Free(vector);
            PyErr_SetString(PyExc_TypeError, "3D vector element items must be float.");
            Py_DECREF(pItem);
            return NULL;
        }
        PyObject_AsDouble( pItem, &(vector[i]) );
        Py_DECREF(pItem);
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
    for (i = 0 ; i < n ; i++) {
        pItem = PyList_GetItem(pListChainID, i);
        Py_INCREF(pItem);

        objectsRepresentation = PyObject_Repr(pItem); //Now a unicode object, new ref
        Py_DECREF(pItem);
        
        PyObject* pyStr = PyUnicode_AsUTF8String(objectsRepresentation); // NEw ref
        Py_DECREF(objectsRepresentation);
        
        s = PyBytes_AS_STRING(pyStr);
        Py_DECREF(pyStr);

        (*buffer)[i] = s[1];

    }
    
    return 1;
}

int unpackString(PyObject *pListOfStrings, char ***buffer) {
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
        Py_INCREF(pItem);

        objectsRepresentation = PyObject_Repr(pItem); //Now a unicode object
        Py_DECREF(pItem);

        PyObject* pyStr = PyUnicode_AsUTF8String(objectsRepresentation);
        Py_DECREF(objectsRepresentation);

        s = PyBytes_AS_STRING(pyStr);
        Py_DECREF(pyStr);

        sLen =  strlen(s); // This corresponds to the actual string surrounded by \' , ie : 'MYTSRING'
       
        (*buffer)[i] = PyMem_New(char, sLen - 1);
        for (int j = 1 ; j < sLen - 1 ; j++) {
            (*buffer)[i][j - 1] = s[j];
        }
        (*buffer)[i][sLen - 2] = '\0';
    }

    return 1;
}

double *unpackCoordinates(PyObject *pListCoor) {
    PyObject *pItem;
    Py_ssize_t n = PyList_Size(pListCoor);
    double *buffer = PyMem_New(double, n);

    for (int i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListCoor, i);
        Py_INCREF(pItem);
        
        if(!PyFloat_Check(pItem)) {
            PyErr_SetString(PyExc_TypeError, "coordinate items must be float.");
            PyMem_Free(buffer);
            Py_DECREF(pItem); 
            return NULL;
        }

        PyObject_AsDouble(pItem, &(buffer[i]) );
        Py_DECREF(pItem);        
    }
    return buffer;
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


/*
    WE NEED ERROR MANAGMENT
*/
atom_t *structDictToAtoms(PyObject *pyDictObject, int *nAtoms) {

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

    /* Create data structures and compute */
    atom_t *atomList = readFromArrays(*nAtoms, coorX, coorY, coorZ, chainID, resSeq, resName, atomName);

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
    fprintf(stderr, "Exiting from structDictToAtoms\n");
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

#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
#include "mesh_io.h"
#include "transform_mesh.h"
#include "encode.h"

struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

//static struct module_state _state;

static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
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

static int unpackVector3(PyObject *pyObject_Tuple, double (*vector)[3]) {
#ifdef DEBUG
    PySys_WriteStdout("--->Unpack Vector3\n");
#endif

    float x;
    PyObject *pItem;
    Py_ssize_t n;
    n = PyTuple_Size(pyObject_Tuple);
    for (int i = 0 ; i < 3 ; i++) {
        pItem = PyTuple_GetItem(pyObject_Tuple, i);
        if(!PyFloat_Check(pItem)) {
            PyErr_SetString(PyExc_TypeError, "3D vector element items must be float.");
            return 0;
        }
        PyObject_AsDouble( pItem, &(*vector)[i]);
    }
    return 1;
}

static int unpackChainID(PyObject *pListChainID, char **buffer) {
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

static int unpackString(PyObject *pListOfStrings, char ***buffer) {
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

static int unpackCoordinates(PyObject *pListCoor, double **buffer) {
#ifdef DEBUG
    PySys_WriteStdout("*** Unpack coordinates ***\n");
#endif
    PyObject *pItem;
    Py_ssize_t n;

    int i;
    n = PyList_Size(pListCoor);
#ifdef DEBUG
    PySys_WriteStdout("Allocating for %d coordinates\n", (int)n);
#endif
    *buffer = PyMem_New(double, n);

    double u;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListCoor, i);
        if(!PyFloat_Check(pItem)) {
            PyErr_SetString(PyExc_TypeError, "coordinate items must be float.");
            PyMem_Free(*buffer);
            return 0;
        }

        PyObject_AsDouble(pItem, &u);
        PyObject_AsDouble(pItem, &(*buffer)[i]);
   //     PySys_WriteStdout("TEST:: %.2f\n", (*buffer)[i] );
    }
    #ifdef DEBUG
    PySys_WriteStderr("Allocation done\n");
    #endif
    return 1;
}

static void freeBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n) {
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

static atom_t *structDictToAtoms(PyObject *pyDictObject, int *nAtoms) {
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
    double *coorX, *coorY, *coorZ;
    unpackCoordinates(pyObj_x, &coorX);
    unpackCoordinates(pyObj_y, &coorY);
    unpackCoordinates(pyObj_z, &coorZ);
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

// Passing a list of two dictionaries, the Euler triplet and the translation triplet
static PyObject *ccmap_compute_zdock_pose(PyObject *self, PyObject *args) {

    PyObject *pyMolStrucTuple, *eulerTuple, *translationTuple, *recOffset, *ligOffset;
    float userThreshold;
    double eulerAngle[3]  = { 0.0, 0.0, 0.0 }; // euler angle to generate pose
    double translation[3] = { 0.0, 0.0, 0.0 }; // translation vector to generate pose
    double offsetLIG[3]   = { 0.0, 0.0, 0.0 }; // translation to center ligand
    double offsetREC[3]   = { 0.0, 0.0, 0.0 }; // translation to center receptor

    int nAtomsRec, nAtomsLig;
    int *ccmap = NULL;
    unsigned int finalLen=0;
    atom_t *atomListRec, *atomListLig;

    if (!PyArg_ParseTuple(args, "O!fO!O!O!O!", &PyTuple_Type, &pyMolStrucTuple, &userThreshold,
                                               &PyTuple_Type, &eulerTuple, &PyTuple_Type, &translationTuple,
                                               &PyTuple_Type, &recOffset, &PyTuple_Type, &ligOffset)) {
    //if (!PyArg_ParseTuple(args, "fo!", &PyList_Type, &userThreshold, &pyDictList)) {
        PyErr_SetString(PyExc_TypeError, "Improper set of arguments °\n");
        return NULL;
    }

    PyObject *pStructAsDictRec = PyTuple_GetItem(pyMolStrucTuple, 0);
    PyObject *pStructAsDictLig = PyTuple_GetItem(pyMolStrucTuple, 1);
    atomListRec = structDictToAtoms(pStructAsDictRec, &nAtomsRec);
    atomListLig = structDictToAtoms(pStructAsDictLig, &nAtomsLig);

    if ( ! unpackVector3(eulerTuple, &eulerAngle) ) {
         PyErr_SetString(PyExc_TypeError, "Fail to unpack euler triplet");
         return NULL;
    }
    if ( ! unpackVector3(translationTuple, &translation) ) {
         PyErr_SetString(PyExc_TypeError, "Fail to unpack translation triplet");
         return NULL;
    }
    if ( ! unpackVector3(ligOffset, &offsetLIG) ) {
         PyErr_SetString(PyExc_TypeError, "Fail to unpack ligand offset triplet");
         return NULL;
    }
    if ( ! unpackVector3(recOffset, &offsetREC) ) {
         PyErr_SetString(PyExc_TypeError, "Fail to unpack receptor offset triplet");
         return NULL;
    }
#ifdef DEBUG
    PySys_WriteStdout("Unpacking euler %.2f %.2f %.2f||translation vectors %.2f %.2f %.2f\n", \
        eulerAngle[0], eulerAngle[1], eulerAngle[2], \
        translation[0], translation[1], translation[2]);
    //  PySys_WriteStdout("Unpacking %d rotation matrices [contact distance is %f]\n", (int)nRotMatrix, userThreshold);
#endif

    Py_BEGIN_ALLOW_THREADS
    transformAtomList(atomListRec, NULL, offsetREC);
    transformAtomList(atomListLig, eulerAngle, offsetLIG);
    transformAtomList(atomListLig, NULL, translation);

    // TO REPLACE HERE
    //ccmap = residueContactMap_DUAL(atomListRec, nAtomsRec, atomListLig, nAtomsLig, userThreshold, &finalLen);
    // printTable(ccmap, finalLen);
    Py_END_ALLOW_THREADS
    PyObject* python_ccmap= PyList_New(finalLen);
    // printf("Size of list %ld \n", PyList_GET_SIZE(python_ccmap));
    for (Py_ssize_t i=0; i<finalLen; i++){
      PyObject *py_value= NULL;
      py_value = Py_BuildValue("i",ccmap[i]);
      int res=PyList_SetItem(python_ccmap, i, py_value );
    }
  
    if ( ! backMapCoordinates(atomListRec, pStructAsDictRec) ) {
        PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
        return NULL;
    }
    if ( ! backMapCoordinates(atomListLig, pStructAsDictLig) ) {
        PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
        return NULL;
    }

    destroyAtomList(atomListRec, nAtomsRec);
    destroyAtomList(atomListLig, nAtomsLig);
    free(ccmap);

    return python_ccmap;
}

static PyObject *
ccmap_compute_ext(PyObject *self, PyObject *args)
{
    char dummy[] = "XY\0";

/*
PyArg_ParseTuple() returns true (nonzero) if all arguments have the right type and its components have been stored in
the variables whose addresses are passed. It returns false (zero) if an invalid argument list was passed.
In the latter case it also raises an appropriate exception so the calling function can return NULL immediately (as we saw in the example).
*/

PyObject *pyDictList, *pTuple;
float userThreshold;
if (!PyArg_ParseTuple(args, "O!f", &PyList_Type, &pyDictList, &userThreshold)) {
    PyErr_SetString(PyExc_TypeError, "parameters must be a list of dictionnaries and a distance value.");
    return NULL;
}
    Py_ssize_t nStructPairs = PyList_Size(pyDictList);
#ifdef DEBUG
    PySys_WriteStdout("Unpacking %d pairs of structure coordinates [contact distance is %f]\n", (int)nStructPairs, userThreshold);
#endif


    //A list to store results

    /*int PyList_SetItem(PyObject *list, Py_ssize_t index, PyObject *item)
    Set the item at index index in list to item. Return 0 on success or -1 on failure.

    Note This function “steals” a reference to item and discards a reference to an item already in the list at the affected position.
    */
    PyObject *PyList_results =  PyList_New(nStructPairs);
    PyObject *pStructAsDictRec, *pStructAsDictLig;

    atom_t *atomListRec, *atomListLig;
    int nAtomsRec, nAtomsLig;
    int *ccmap = NULL;
    unsigned int finalLen=0;

    //return Py_BuildValue("s", dummy);


    for (int i = 0; i < (int)nStructPairs ; i++) {
#ifdef DEBUG
        PySys_WriteStdout("____Structure_____ %d\n", i);
#endif

        pTuple = PyList_GetItem(pyDictList, i);
        pStructAsDictRec = PyTuple_GetItem(pTuple, 0);
        pStructAsDictLig = PyTuple_GetItem(pTuple, 1);

#ifdef DEBUG
        PySys_WriteStdout("REF COUNT pStructAsDictRec :: is %d\nREF COUNT pStructAsDictLig :: is %d\n", Py_REFCNT(pStructAsDictRec), Py_REFCNT(pStructAsDictLig));
#endif

        atomListRec = structDictToAtoms(pStructAsDictRec, &nAtomsRec);
        atomListLig = structDictToAtoms(pStructAsDictLig, &nAtomsLig);

#ifdef DEBUG
        PySys_WriteStdout("USED__REF COUNT pStructAsDictRec :: is %d\nUSED__REF COUNT pStructAsDictLig :: is %d\n", Py_REFCNT(pStructAsDictRec), Py_REFCNT(pStructAsDictLig));
#endif
        // TO REPLACE HERE
        //ccmap = residueContactMap_DUAL(atomListRec, nAtomsRec, atomListLig, nAtomsLig, userThreshold, &finalLen);
        PyObject* python_ccmap= PyList_New(finalLen);
    // printf("Size of list %ld \n", PyList_GET_SIZE(python_ccmap));
        for (Py_ssize_t i=0; i<finalLen; i++){
      // printf(" %p  ; %ld ; %d ", python_ccmap, i, ccmap[i] );
            PyObject *py_value= NULL;
            py_value = Py_BuildValue("i",ccmap[i]);
            int res=PyList_SetItem(python_ccmap, i, py_value );
      // printf("%d \n", res);
        }

        PyList_SetItem(PyList_results, i, Py_BuildValue("O", python_ccmap));

        //PyList_SetItem(PyList_results, i, Py_BuildValue("s", "TOTOTO"));
#ifdef DEBUG
        PySys_WriteStderr("Destroying Atoms lists\n");
#endif

        destroyAtomList(atomListRec, nAtomsRec);
        destroyAtomList(atomListLig, nAtomsLig);

#ifdef DEBUG
        PySys_WriteStdout("ATOM_DESTROYED__REF COUNT pStructAsDictRec :: is %d\nATOM_DESTROYED__REF COUNT pStructAsDictLig :: is %d\n", Py_REFCNT(pStructAsDictRec), Py_REFCNT(pStructAsDictLig));
        PySys_WriteStdout("Freeing json C pointer safely\n");
#endif

        free(ccmap);

#ifdef DEBUG
        PySys_WriteStdout("CCMAP_AND_ATOM_DESTROYED__REF COUNT pStructAsDictRec :: is %d\nCCMAP_AND_ATOM_DESTROYED__REF COUNT pStructAsDictLig :: is %d\n", Py_REFCNT(pStructAsDictRec), Py_REFCNT(pStructAsDictLig));
#endif

        // HACK -- makes program work //
        // Delalocation of the passed arguments makes
        // decrements the pStructAsDict??? refernces count
        //Py_INCREF(pStructAsDictRec);
        //Py_INCREF(pStructAsDictLig);
    }
#ifdef DEBUG
    PySys_WriteStderr("Going out\n");
    PySys_WriteStdout("REF COUNT results  :: is %d\n", Py_REFCNT(PyList_results) );
    PySys_WriteStdout("REF COUNT pDict  :: is %d\n", Py_REFCNT(pyDictList) );
    PySys_WriteStdout("REF COUNT pDictRec :: is %d\n", Py_REFCNT(pStructAsDictRec) );
    PySys_WriteStdout("REF COUNT pDictLig  :: is %d\n", Py_REFCNT(pStructAsDictLig) );
#endif

   return PyList_results;
   // return Py_BuildValue("s", dummy);
}

static PyObject *ccmap_compute_base(PyObject *self, PyObject *args) {

PyObject *pListX, *pListY, *pListZ, *pListChainID, *pListResName, *pListResSeq, *pListAtomName, *pComputingType, *atomicBool;
atomicBool = NULL;
float userThreshold;
if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!f|O", &PyList_Type, &pListX, \
                                                 &PyList_Type, &pListY, \
                                                 &PyList_Type, &pListZ, \
                                                 &PyList_Type, &pListChainID,\
                                                 &PyList_Type, &pListResSeq, \
                                                 &PyList_Type, &pListResName, \
                                                 &PyList_Type, &pListAtomName, \
                                                 &userThreshold, \
                                                 &atomicBool )) {
    PyErr_SetString(PyExc_TypeError, "Parameters must be Seven and a distance value lists.");
    return NULL;
}
    bool bAtomic = false;
    if (atomicBool == NULL) {
        //PySys_WriteStdout("##applying default bAtom");
    } else {
        //PySys_WriteStdout("processing bAtom");
        if (PyObject_IsTrue(atomicBool))
            bAtomic = true;
        Py_XDECREF(atomicBool);
    }

    Py_ssize_t n = PyList_Size(pListX);
#ifdef DEBUG
    PySys_WriteStdout("Unpacking %d atoms [contact distance is %f]\n", (int)n, userThreshold);
#endif
    /* Unpacking atom specification vectors */
    double *coorX, *coorY, *coorZ;
    unpackCoordinates(pListX, &coorX);
    unpackCoordinates(pListY, &coorY);
    unpackCoordinates(pListZ, &coorZ);

    char *chainID;
    unpackChainID(pListChainID, &chainID);

    char **resSeq;
    unpackString(pListResSeq, &resSeq);

    char **resName;
    unpackString(pListResName, &resName);

    char **atomName;
    unpackString(pListAtomName, &atomName);
    /* Create data structures and compute */
    atom_t *atomList = readFromArrays(n, coorX, coorY, coorZ, chainID, resSeq, resName, atomName);
    // Old API
    //char *ccmap = residueContactMap(atomList, n, userThreshold);
    // New API
    bool bEncode = false;
   // bool bAtomic = false;
    ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                        ? &atomicContactMap
                        : &residueContactMap;
    ccmapView_t *ccmapView = computeMap(atomList, n, NULL, 0, userThreshold, bEncode);
    
    PyObject *rValue = Py_BuildValue("s", ccmapView->asJSON);
    destroyCcmapView(ccmapView);
    
    return rValue;
}

/*
static char* keywords[] = {"", "", "k1", "k2", "k3", "k4", NULL};

// [...]

PyArg_ParseTupleAndKeywords(args, kwargs,
                            "O!|O!$O!idO", keywords,
                            &PyArray_Type, &arg1, &PyArray_Type, &arg2,
                            &PyArray_Type, &k1, &k2, &k3, &k4);
*/
// ccmap(coordinatesDictA, y=coordinatesDictB, value=4.5, atomic=False)
// Print refdict count on entry / exit
static PyObject *ccmap_compute_flex(PyObject* self, PyObject* args, PyObject* kwargs) {
    PySys_WriteStderr("ccmap_compute_flex");
    static char *kwlist[] = {"", "y", "dist", "atomic", "encode", NULL};
/*PyObject *pListX_i, *pListY_i, *pListZ_i, *pListChainID_i, *pListResName_i, *pListResSeq_i, *pListAtomName_i;
PyObject *pListX_j, *pListY_j, *pListZ_j, *pListChainID_j, *pListResName_j, *pListResSeq_j, *pListAtomName_j;

 PyObject *pComputingType, *atomicBool;
*/
PyObject *atomicBool = NULL;
PyObject *encodeBool = NULL;
float userThreshold = 4.5;
PyObject *coorDictI;
PyObject *coorDictJ = NULL;
if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
                                "O!|O!fOO", kwlist,\
                                &PyDict_Type, &coorDictI, \
                                &PyDict_Type, &coorDictJ, \
                                &userThreshold, \
                                &atomicBool,
                                &encodeBool )) {
    PySys_WriteStderr("#####");
    PyErr_SetString(PyExc_TypeError, "Parameters must be One dict and one dict optional distance value and bAtomic.");
    return NULL;
}

// Parsing map type and encoding scheme
bool bAtomic = false;
bool bEncode = false;
if (atomicBool != NULL) {
    if (PyObject_IsTrue(atomicBool))
        bAtomic = true;
    Py_XDECREF(atomicBool);
}
if (encodeBool != NULL) {
    if (PyObject_IsTrue(encodeBool))
        bEncode = true;
    Py_XDECREF(encodeBool);
}

// Parsing coordinates
int iLen, jLen;
atom_t *iAtomList = structDictToAtoms(coorDictI, &iLen);
atom_t *jAtomList = NULL;
if(coorDictJ != NULL)
    jAtomList = structDictToAtoms(coorDictJ, &jLen);
PySys_WriteStderr("Running with %s coordinates sets dist=%f bAtom:%s bEncode:%s",\
                                    jAtomList == NULL ? "one" : "two",\
                                    userThreshold,\
                                    bAtomic ? "true" : "false",\
                                    bEncode ? "true" : "false");

// Computing
ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                    ? &atomicContactMap
                    : &residueContactMap;
ccmapView_t *ccmapView = computeMap(iAtomList, iLen, jAtomList, jLen, userThreshold, bEncode);

// Build returned pyObject
PyObject *rValue;
if (bEncode) { // ccmap is int Vec<int> form, we build/return an pyList from it
    rValue = PyList_New(ccmapView->encodeLen);
    for (Py_ssize_t i = 0; i < ccmapView->encodeLen; i++){
        PyObject *py_value = NULL;
        py_value = Py_BuildValue("i", ccmapView->asENCODE[i]);
        PyList_SetItem(rValue, i, py_value);
    }
       // PyList_SetItem(PyList_results, i, Py_BuildValue("O", python_ccmap));
    } else {
    rValue = Py_BuildValue("s", ccmapView->asJSON);
}
// Cleaning 
destroyAtomList(iAtomList, iLen);
if(coorDictJ != NULL)
    destroyAtomList(jAtomList, jLen);
destroyCcmapView(ccmapView);


return rValue;

/*
PyObject* xCoor = PyDict_GetItemString(coorDictI, "x");
PyObject* yCoor = PyDict_GetItemString(coorDictI, "y");
PyObject* zCoor = PyDict_GetItemString(coorDictI, "z");

if (xCoor == NULL){ 
    PySys_WriteStderr("no x key found");
    PyErr_SetString(PyExc_TypeError, "no x key found");
    return NULL;
}
if ( !PyList_Check(xCoor) ) {
    PyErr_SetString(PyExc_TypeError, "x value is not a list");
    return NULL;
}
if (yCoor == NULL){ 
    PySys_WriteStderr("no y key found");
    PyErr_SetString(PyExc_TypeError, "no y key found");
    return NULL;
}
if ( !PyList_Check(yCoor) ) {
    PyErr_SetString(PyExc_TypeError, "y value is not a list");
    return NULL;
}
if (zCoor == NULL){ 
    PySys_WriteStderr("no z key found");
    PyErr_SetString(PyExc_TypeError, "no z key found");
    return NULL;
}
if ( !PyList_Check(zCoor) ) {
    PyErr_SetString(PyExc_TypeError, "z value is not a list");
    return NULL;
}

Py_ssize_t coorLen = PyList_Size(xCoor);
if(coorLen != PyList_Size(yCoor) || coorLen != PyList_Size(zCoor)) {
    PyErr_SetString(PyExc_TypeError, "uneven x,y,z coordinates list");
    return NULL;
}

double *coorX, *coorY, *coorZ;
unpackCoordinates(xCoor, &coorX);
unpackCoordinates(yCoor, &coorY);
unpackCoordinates(zCoor, &coorZ);

for(int i = 0; i < coorLen; i++) {
    PySys_WriteStdout("%g %g %g\n", coorX[i], coorY[i], coorZ[i]);
}
*/

/*

    char dummy[] = "Dval";

     PySys_WriteStdout("value=%g", userThreshold);
return Py_BuildValue("s", dummy);
*/
/*
    bool bAtomic = false;
    if (atomicBool == NULL) {
        //PySys_WriteStdout("##applying default bAtom");
    } else {
        //PySys_WriteStdout("processing bAtom");
        if (PyObject_IsTrue(atomicBool))
            bAtomic = true;
        Py_XDECREF(atomicBool);
    }

    if (coorDictJ == NULL) {
        PySys_WriteStdout("No second dict");
    } else{
        PySys_WriteStdout("Found second dict");
    }
    */
    /*
    Py_ssize_t n = PyList_Size(pListX);
#ifdef DEBUG
    PySys_WriteStdout("Unpacking %d atoms [contact distance is %f]\n", (int)n, userThreshold);
#endif
    double *coorX, *coorY, *coorZ;
    unpackCoordinates(pListX, &coorX);
    unpackCoordinates(pListY, &coorY);
    unpackCoordinates(pListZ, &coorZ);

    char *chainID;
    unpackChainID(pListChainID, &chainID);

    char **resSeq;
    unpackString(pListResSeq, &resSeq);

    char **resName;
    unpackString(pListResName, &resName);

    char **atomName;
    unpackString(pListAtomName, &atomName);
    atom_t *atomList = readFromArrays(n, coorX, coorY, coorZ, chainID, resSeq, resName, atomName);
    bool bEncode = false;
    ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                        ? &atomicContactMap
                        : &residueContactMap;
    ccmapView_t *ccmapView = computeMap(atomList, n, NULL, 0, userThreshold, bEncode);
    
    PyObject *rValue = Py_BuildValue("s", ccmapView->asJSON);
    destroyCcmapView(ccmapView);
    
    return rValue;
    */
}


int PyObject_AsDouble(PyObject *py_obj, double *x)
{
  PyObject *py_float;

  py_float = PyNumber_Float(py_obj);

  if (py_float == NULL) return -1;

  *x = PyFloat_AsDouble(py_float);

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

/* MODULE DECLARATION */



static PyMethodDef CcmapMethods[] = {
    //...
    {"compute",  ccmap_compute_base, METH_VARARGS,
     "Compute a residue or atomic contact map."},
     {"dvl",   ccmap_compute_flex, METH_VARARGS | METH_KEYWORDS,
     "Compute a residue or atomic contact map."},
    {"duals",  ccmap_compute_ext, METH_VARARGS,
     "Compute a protein-protein interface residue contact map."},
     {"zmap",  ccmap_compute_zdock_pose, METH_VARARGS,
     "Compute a contact map from a zdock-like encoded pose"},
    //...
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


static int myextension_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int myextension_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "ccmap",
        NULL,
        sizeof(struct module_state),
        CcmapMethods,
        NULL,
        myextension_traverse,
        myextension_clear,
        NULL
};
#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_ccmap(void)
{
    PyObject *module = PyModule_Create(&moduledef);

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("myextension.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    return module;
}


//
//int
//main(int argc, char *argv[])
//{
    /* Pass argv[0] to the Python interpreter */
//    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
//    Py_Initialize();

    /* Add a static module */
//    init_ccmap();
//}

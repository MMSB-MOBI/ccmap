#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
#include "mesh_io.h"
#include "transform_mesh.h"

/*
    Python C-API, pure-C  implementation interface
*/

/*  GL 16 / 10 / 2017
This function is still in progress and disabled in released module

// We use Py_Mem and need to hold the GIL
static
void allocateList(char **ccmapList, atom_t **atomListRecList, atomt_t **atomListRecList, double **eulerAngleList, double **translationList, Py_ssize_t nComplexes) {
    int n = int(nComplexes);
    ccmapList = PyMem_New(*char, nComplexes);
    atomListRecList = PyMem_New(*atom_t, nComplexes);
    atomListLigList = PyMem_New(*atom_t, nComplexes);
    eulerAngleList  = PyMem_New(*double, nComplexes);
    translationList = = PyMem_New(*double, nComplexes);
    for (int i = 0; i < n ; i++) {
        eulerAngleList[i] = PyMem_New(*double, 3);
        translationList[i] = PyMem_New(*double, 3);
    }
}

static freeList(char **ccmapList, atom_t **atomListRecList, atomt_t **atomListRecList, Py_ssize_t  nComplexes,
                int nRecAtoms, int nLigAtoms) {
    int n = int(nComplexes);
    for(int i = 0; i < n; i++) {
        destroyAtomList(atomListRecList[i], nRecAtoms);
        destroyAtomList(atomListLigList[i], nLigAtoms);
        free(ccmapList[i]);
        PyMem_Free(eulerAngleList[i]);
        PyMem_Free(TranslationList[i]);
    }
    PyMem_Free(ccmapList);
    PyMem_Free(atomListRecList);
    PyMem_Free(atomListLigList);
    PyMem_Free(eulerAngleList);
    PyMem_Free(translationList);
}
*/

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
            //PySys_WriteStdout("TOTO\n");
            //pItem = PyMem_New(double, 1);

            //PySys_WriteStdout("TOTO :: %f<<\n", atomCurrent->x);


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
   // PySys_WriteStdout("--->%d\n", n);

    *buffer = PyMem_New(char, n);

    /*PyObject* objectsRepresentation = PyObject_Repr(yourObject);
    const char* s = PyString_AsString(objectsRepresentation);
    */
    PyObject* objectsRepresentation;
    const char* s;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListChainID, i);
        objectsRepresentation = PyObject_Repr(pItem);
        s = PyString_AsString(objectsRepresentation);
        (*buffer)[i] = s[1];
    }
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

    /*PyObject* objectsRepresentation = PyObject_Repr(yourObject);
    const char* s = PyString_AsString(objectsRepresentation);
    */
    PyObject* objectsRepresentation;
    const char* s;
    int sLen;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListOfStrings, i);
        objectsRepresentation = PyObject_Repr(pItem);
        s = PyString_AsString(objectsRepresentation); // DOC says it must not be de allocated
        sLen =  strlen(s); // This corresponds to the actual string surrounded by \' , ie : 'MYTSRING'
        //PySys_WriteStdout("--->%s[%d]\n", s, strlen(s));
        (*buffer)[i] = PyMem_New(char, sLen - 1);
        for (int j = 1 ; j < sLen - 1 ; j++) {
            (*buffer)[i][j - 1] = s[j];
        }
        (*buffer)[i][sLen - 2] = '\0';
        Py_DECREF(objectsRepresentation);
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
static PyObject *
ccmap_compute_zdock_pose(PyObject *self, PyObject *args) {

    PyObject *pyMolStrucTuple, *eulerTuple, *translationTuple, *recOffset, *ligOffset;
    float userThreshold;
    double eulerAngle[3]  = { 0.0, 0.0, 0.0 }; // euler angle to generate pose
    double translation[3] = { 0.0, 0.0, 0.0 }; // translation vector to generate pose
    double offsetLIG[3]   = { 0.0, 0.0, 0.0 }; // translation to center ligand
    double offsetREC[3]   = { 0.0, 0.0, 0.0 }; // translation to center receptor

    int nAtomsRec, nAtomsLig;
    char *ccmap = NULL;
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


    ccmap = residueContactMap_DUAL(atomListRec, nAtomsRec, atomListLig, nAtomsLig, userThreshold);
    Py_END_ALLOW_THREADS
   // PyList_SetItem(PyList_results, i, Py_BuildValue("s", ccmap));

    PyObject *pyString_ccmap = Py_BuildValue("s", ccmap);

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



    return pyString_ccmap;
}

// To fill up the thread we develop a function similar to ccmap_compute_zdock_pose
// Which performs several transforamtion in a row
// Passing a list of two dictionaries, AND a LIST OF Euler triplets and a LIST OF the translation triplet

/*  GL 16 / 10 / 2017
This function is still in progress and disabled in released module

static PyObject *
ccmap_compute_zdock_poseList(PyObject *self, PyObject *args) {

    PyObject *pyMolStrucTuple, *eulerTupleList, *translationTupleList, *recOffset, *ligOffset, *pTuple;
    float userThreshold;
    double **eulerAngleList; // list of euler angle to generate pose
    double **translationList; // translation vector to generate pose
    double offsetLIG[3]   = { 0.0, 0.0, 0.0 }; // translation to center ligand
    double offsetREC[3]   = { 0.0, 0.0, 0.0 }; // translation to center receptor

    int nAtomsRec, nAtomsLig;
    char **ccmapList = NULL;
    atom_t *atomListRec, *atomListLig;

    if (!PyArg_ParseTuple(args, "O!fO!O!O!O!", &PyList_Type, &pyMolStrucTupleList, &userThreshold,
                                               &PyList_Type, &eulerTupleList, &PyList_Type, &translationTupleList,
                                               &PyTuple_Type, &recOffset, &PyTuple_Type, &ligOffset)) {
    //if (!PyArg_ParseTuple(args, "fo!", &PyList_Type, &userThreshold, &pyDictList)) {
        PyErr_SetString(PyExc_TypeError, "Improper set of arguments °\n");
        return NULL;
    }

    if ( ! unpackVector3(ligOffset, &offsetLIG) ) {
         PyErr_SetString(PyExc_TypeError, "Fail to unpack ligand offset triplet");
         return NULL;
    }
    if ( ! unpackVector3(recOffset, &offsetREC) ) {
         PyErr_SetString(PyExc_TypeError, "Fail to unpack receptor offset triplet");
         return NULL;



    Py_ssize_t nStructPairs = PyList_Size(pyMolStrucTupleList);

    // malloc atomListRecList
    // malloc atomListLigList

    allocateList(ccmapList, atomListRecList, atomListRecList, eulerAngleList, translationList, nStructPairs);
    PyObject *PyList_results =  PyList_New(nStructPairs);

    for (int i = 0; i < (int)nStructPairs ; i++) {

        pTuple = PyList_GetItem(pyDictList, i);
        pStructAsDictRec = PyTuple_GetItem(pTuple, 0);
        pStructAsDictLig = PyTuple_GetItem(pTuple, 1);
        PyObject *pStructAsDictRec = PyTuple_GetItem(pyMolStrucTuple, 0);
        PyObject *pStructAsDictLig = PyTuple_GetItem(pyMolStrucTuple, 1);
        atomListRecList[i] = structDictToAtoms(pStructAsDictRec, &nAtomsRec);
        atomListLigList[i] = structDictToAtoms(pStructAsDictLig, &nAtomsLig);
        if ( ! unpackVector3(eulerTuple, &eulerAngleList[i]) ) {
             PyErr_SetString(PyExc_TypeError, "Fail to unpack euler triplet");
            return NULL;
        }
        if ( ! unpackVector3(translationTuple, &translationList[i]) ) {
             PyErr_SetString(PyExc_TypeError, "Fail to unpack translation triplet");
            return NULL;
        }

#ifdef DEBUG
        PySys_WriteStdout("Unpacking euler %.2f %.2f %.2f||translation vectors %.2f %.2f %.2f\n", \
            eulerAngle[0], eulerAngle[1], eulerAngle[2], \
            translation[0], translation[1], translation[2]);
    //  PySys_WriteStdout("Unpacking %d rotation matrices [contact distance is %f]\n", (int)nRotMatrix, userThreshold);
#endif
    }

    Py_BEGIN_ALLOW_THREADS
    for (int i = 0; i < (int)nStructPairs ; i++) {
        transformAtomList(atomListRecList[i], NULL, offsetREC);
        transformAtomList(atomListLigList[i], eulerAngleList[i], offsetLIG);
        transformAtomList(atomListLigList[i], NULL, translationList[i]);
        ccmapList[i] = residueContactMap_DUAL(atomListRecList[i], nAtomsRec, atomListLigList[i], nAtomsLig, userThreshold);
    }
    Py_END_ALLOW_THREADS

    for (int i = 0; i < (int)nStructPairs ; i++) {
        if ( ! backMapCoordinates(atomListRec, pStructAsDictRec) ) {
            PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
            return NULL;
        }
        if ( ! backMapCoordinates(atomListLig, pStructAsDictLig) ) {
            PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
            return NULL;
        }
        PyList_SetItem(PyList_results, i, Py_BuildValue("s", ccmapList[i]));
    }


    freeList(ccmapList, atomListRecList, atomListRecList, eulerAngleList, translationList, nStructPairs, nAtomsRec, nAtomsLig);

    return PyList_results;
}
*/
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
    char *ccmap = NULL;


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

        ccmap = residueContactMap_DUAL(atomListRec, nAtomsRec, atomListLig, nAtomsLig, userThreshold);
        PyList_SetItem(PyList_results, i, Py_BuildValue("s", ccmap));

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


static PyObject *ccmap_compute(PyObject *self, PyObject *args) {
    char dummy[] = "XY\0";

/*
PyArg_ParseTuple() returns true (nonzero) if all arguments have the right type and its components have been stored in
the variables whose addresses are passed. It returns false (zero) if an invalid argument list was passed.
In the latter case it also raises an appropriate exception so the calling function can return NULL immediately (as we saw in the example).
*/

PyObject *pListX, *pListY, *pListZ, *pListChainID, *pListResName, *pListResSeq, *pListAtomName;
float userThreshold;
if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!f", &PyList_Type, &pListX, &PyList_Type, &pListY, &PyList_Type, &pListZ, &PyList_Type, &pListChainID, &PyList_Type, &pListResSeq, &PyList_Type, &pListResName, &PyList_Type, &pListAtomName,&userThreshold)) {
    PyErr_SetString(PyExc_TypeError, "parameters must be Seven and a distance value lists.");
    return NULL;
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

    char *ccmap = residueContactMap(atomList, n, userThreshold);

    atomList = destroyAtomList(atomList, n);

    return Py_BuildValue("s", ccmap);
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

static PyMethodDef CcmapMethods[] = {
    //...
    {"compute",  ccmap_compute, METH_VARARGS,
     "Compute a contact map."},
    {"duals",  ccmap_compute_ext, METH_VARARGS,
     "Scores interface based on their contact map."},
     {"zmap",  ccmap_compute_zdock_pose, METH_VARARGS,
     "Compute a contact map from a zdock-like encoded pose"},
    //...
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initccmap(void)
{
    (void) Py_InitModule("ccmap", CcmapMethods);
}


int
main(int argc, char *argv[])
{
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initccmap();
}

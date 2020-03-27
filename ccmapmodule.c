#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
#include "mesh_io.h"
#include "transform_mesh.h"
#include "encode.h"
#include "ccmapmodule_utils.h"



PyDoc_STRVAR(
    ccmap_compute_zdock_pose_doc,
    "foo(timeout, flags=None, /)\n"
    "--\n"
    "\n"
    "Great example function\n"
    "Arguments: (timeout, flags=None)\n"
    "Doc blahblah doc doc doc.");



struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
/*
static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}
*/
/*
Python API : "zmap"
parameters : 2-tuple of dictorizedStructures
             contact distance float
             3-uple of euler angles, aka docking pose-specific rotation
             3-uple of translation of one translation vector components, aka docking-pose specific translation
             3-uple of translation of one translation vector components, aka receptor translation offset
             3-uple of translation of one translation vector components, aka ligand translation offset
             Encoding flag, Boolean
Returns: a string if Encoding flag is false, a list of residue numbers otherwise
*/

// Passing a list of two dictionaries, the Euler triplet and the translation triplet
static PyObject *ccmap_compute_zdock_pose(PyObject *self, PyObject *args) 
{
    PySys_WriteStderr("COUCOU\n");
PyObject *pyMolStrucTuple, *eulerTuple, *translationTuple, *recOffset, *ligOffset;
float userThreshold;

int nAtomsRec, nAtomsLig;
atom_t *atomListRec, *atomListLig;
PyObject *encodeBool = NULL;
if (!PyArg_ParseTuple(args, "O!fO!O!O!O!|O", &PyTuple_Type, &pyMolStrucTuple, &userThreshold,
                                            &PyTuple_Type, &eulerTuple, &PyTuple_Type, &translationTuple,
                                            &PyTuple_Type, &recOffset, &PyTuple_Type, &ligOffset, &encodeBool )) {
    PyErr_SetString(PyExc_TypeError, "Improper set of arguments Zmap\n");
    return NULL;
}
bool bEncode = false;
if (encodeBool != NULL) {
    if (PyObject_IsTrue(encodeBool))
        bEncode = true;
    Py_XDECREF(encodeBool);
}
PyObject *pStructAsDictRec = PyTuple_GetItem(pyMolStrucTuple, 0);
PyObject *pStructAsDictLig = PyTuple_GetItem(pyMolStrucTuple, 1);
atomListRec = structDictToAtoms(pStructAsDictRec, &nAtomsRec);
atomListLig = structDictToAtoms(pStructAsDictLig, &nAtomsLig);

double *eulerAngle  = unpackVector3(eulerTuple); // euler angle to generate pose
double *translation = unpackVector3(translationTuple) ; // translation vector to generate pose
double *offsetLIG   = unpackVector3(ligOffset); // translation to center ligand
double *offsetREC   = unpackVector3(recOffset); // translation to center receptor
if (eulerAngle == NULL)
        PyErr_SetString(PyExc_TypeError, "Fail to unpack euler triplet");
if (translation == NULL) 
        PyErr_SetString(PyExc_TypeError, "Fail to unpack translation triplet");        
if (offsetLIG == NULL) 
        PyErr_SetString(PyExc_TypeError, "Fail to unpack ligand offset triplet");
if (recOffset == NULL)
        PyErr_SetString(PyExc_TypeError, "Fail to unpack receptor offset triplet");
if (recOffset == NULL || offsetLIG == NULL || translation == NULL || eulerAngle == NULL) {
    destroyVector3(eulerAngle);
    destroyVector3(translation);
    destroyVector3(offsetLIG);
    destroyVector3(offsetREC);
}    

#ifdef DEBUG
PySys_WriteStdout("Unpacking euler %.2f %.2f %.2f||translation vectors %.2f %.2f %.2f\n", \
    eulerAngle[0], eulerAngle[1], eulerAngle[2], \
    translation[0], translation[1], translation[2]);
//  PySys_WriteStdout("Unpacking %d rotation matrices [contact distance is %f]\n", (int)nRotMatrix, userThreshold);
#endif

ccmapView_t *ccmapView;
Py_BEGIN_ALLOW_THREADS
char titi[200];
stringifyAtom(&atomListRec[0], titi);
fprintf(stderr, "Z1:%s\n", titi);

transformAtomList(atomListRec, NULL, offsetREC);

stringifyAtom(&atomListRec[0], titi);
fprintf(stderr, "Z2:%s\n", titi);

transformAtomList(atomListLig, eulerAngle, offsetLIG);
transformAtomList(atomListLig, NULL, translation);

// Computing
bool bAtomic = false;
ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                ? &atomicContactMap
                : &residueContactMap;
ccmapView = computeMap(atomListRec, nAtomsRec, atomListLig, nAtomsLig, userThreshold, bEncode);
Py_END_ALLOW_THREADS

PyObject *rValue = ccmapViewToPyObject(ccmapView, bEncode);

if ( ! backMapCoordinates(atomListRec, pStructAsDictRec) ) {
    PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
    return NULL;
}
if ( ! backMapCoordinates(atomListLig, pStructAsDictLig) ) {
    PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
    return NULL;
}

destroyVector3(eulerAngle);
destroyVector3(translation);
destroyVector3(offsetLIG);
destroyVector3(offsetREC);
destroyAtomList(atomListRec, nAtomsRec);
destroyAtomList(atomListLig, nAtomsLig);
destroyCcmapView(ccmapView);

return rValue;
}


/*

Python API lzmap

*/

static PyObject *ccmap_compute_zdock_pose_list(PyObject *self, PyObject *args, PyObject* kwargs) 
{
#ifdef DEBUG
PySys_WriteStderr("ccmap_compute_zdock_pose_list");
#endif
    static char *kwlist[] = {"", "", "", "", "offsetRec", "offsetLig", "distance", "encode", NULL};
    PyObject *pyDictRec, *pyDictLig;
    PyObject *eulerList, *transList;
    PyObject *offsetRec = NULL;
    PyObject *offsetLig = NULL;
    float userThreshold = 4.5;
    PyObject *encodeBool = NULL;

    double *offsetLigVector = NULL;
    double *offsetRecVector = NULL;
    
if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
                     "O!O!O!O!|O!O!fO", kwlist, \
                    &PyDict_Type, &pyDictRec,   \
                    &PyDict_Type, &pyDictLig,   \
                    &PyList_Type,  &eulerList,   \
                    &PyList_Type,  &transList,   \
                    &PyList_Type,  &offsetRec,   \
                    &PyList_Type,  &offsetLig,   \
                    &userThreshold,
                    &encodeBool                )) {
    PyErr_SetString(PyExc_TypeError, "Wrong parameters");
    return NULL;
}
if( PyList_Size(eulerList) != PyList_Size(transList) ) {
    PyErr_SetString(PyExc_TypeError, "Transformations tuples lists are of unequal length");
    return NULL;
}

bool bEncode = false;
if (encodeBool != NULL) {
    if (PyObject_IsTrue(encodeBool))
        bEncode = true;
        Py_XDECREF(encodeBool);
}

// Read offset translations for ligand receptor
if (offsetLig != NULL) {  
    offsetLigVector =  unpackVector3(offsetLig);    
    if (offsetLigVector == NULL ) {
        PySys_WriteStderr("Fail to unpack ligand offset triplet");
        return NULL;
    }
}
if (offsetRec != NULL) {  
    offsetRecVector =  unpackVector3(offsetRec);    
    if (offsetRecVector == NULL ) {
        destroyVector3(offsetLigVector);
        PySys_WriteStderr("Fail to unpack receptor offset triplet");
        return NULL;
    }
}

Py_ssize_t nPose;
double **eulerTriplets = createListVector3(eulerList, &nPose);
if(eulerTriplets == NULL) {
    PySys_WriteStderr("Fail to unpack eulers tuples list");
    destroyVector3(offsetLigVector);
    destroyVector3(offsetRecVector);
    return NULL;
}
double **transTriplets = createListVector3(transList, &nPose);
if(transTriplets == NULL) {
    PySys_WriteStderr("Fail to unpack tranlsations tuples list");
    destroyVector3(offsetLigVector);
    destroyVector3(offsetRecVector);
    destroyListVector3(eulerTriplets, nPose);
    return NULL;
}
/*
PySys_WriteStdout("Unpacking %lu  euler and translation poses with bEncode:%s and D=%f\n",\
                                            nPose, bEncode ? "true" : "false", userThreshold);
*/
for (int k = 0; k < nPose ; k++) {
    PySys_WriteStdout("Unpacking %f %f %f // %f %f %f\n", \
                        eulerTriplets[k][0], eulerTriplets[k][1], eulerTriplets[k][2],\
                        transTriplets[k][0], transTriplets[k][1], transTriplets[k][2]\
                        );
}

int nAtomsRec = 0;
int nAtomsLig = 0;
// Reading provided ligand, receptor conformations
atom_t *atomListRec = structDictToAtoms(pyDictRec, &nAtomsRec);
atom_t *atomListLig = structDictToAtoms(pyDictLig, &nAtomsLig);
atom_t *atomListLigBuffer = structDictToAtoms(pyDictLig, &nAtomsLig);

// Alloc for results structures
ccmapView_t **ccmapViews = PyMem_New(ccmapView_t *, nPose); // MAybe not safe has access and subsequent malloc occur in thread below
ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = &residueContactMap;

Py_BEGIN_ALLOW_THREADS
char titi[200];
stringifyAtom(&atomListRec[0], titi);
fprintf(stderr, "LZ1: %s\n", titi);

transformAtomList(atomListRec, NULL, offsetRecVector);

stringifyAtom(&atomListRec[0], titi);
fprintf(stderr, "LZ2: %s\n", titi);
//printAtomList(atomListRec, stdout);
/* Euler rotatio appearently have to be applied b4 ligand offset
if(offsetLigVector != NULL)
    transformAtomList(atomListLig, NULL, offsetLigVector);
*/
// Loop through all poses transformations
for (int iPose = 0 ; iPose < nPose ; iPose++) {
    // Transform Ligand
    // Apparently wrong transforamtion sequence
    //transformAtomList(atomListLigBuffer, eulerTriplets[iPose], transTriplets[iPose]);
    //transformAtomList(atomListLigBuffer, NULL, transTriplets[iPose]);
    transformAtomList(atomListLigBuffer, eulerTriplets[iPose], offsetLigVector);
    
    transformAtomList(atomListLigBuffer, NULL                , transTriplets[iPose]);
    ccmapViews[iPose] = computeMap(atomListRec, nAtomsRec, atomListLigBuffer, nAtomsLig, userThreshold, bEncode);
    // Reset ligand coordinates
    applyCoordinates(atomListLig, atomListLigBuffer);
}
Py_END_ALLOW_THREADS

/* Build return value */
PyObject *rValue = ccmapViewsToPyObject(ccmapViews, nPose, bEncode);

/* transformations cleanup */
destroyVector3(offsetLigVector);
destroyVector3(offsetRecVector);
destroyListVector3(eulerTriplets, nPose);
destroyListVector3(transTriplets, nPose);
/* Computations results cleanup*/
for (int i = 0; i < nPose ; i++)
    destroyCcmapView(ccmapViews[i]);
PyMem_Free(ccmapViews);
/* Coordinates cleanup*/
destroyAtomList(atomListRec, nAtomsRec);
destroyAtomList(atomListLig, nAtomsLig);
destroyAtomList(atomListLigBuffer, nAtomsLig);

return rValue;

    /*
    ccmapView_t ccmapView;
//  Py_BEGIN_ALLOW_THREADS
    transformAtomList(atomListRec, NULL, offsetREC);
    transformAtomList(atomListLig, NULL, translation);
    for (int iPose = 0 ; iPose < )
*/
//  Py_END_ALLOW_THREADS

    //destroyListVector3(eulerTriplets, i);
//    destroyListVector3(transTriplets, j);
    
  //  PyObject *rValue = Py_BuildValue("s", "RTOTO");
   // return rValue;
    
}
/*

    Python API lmap

*/
static PyObject *ccmap_compute_list(PyObject *self, PyObject *args)
{

PyObject *pyDictList, *encodeBool;
float userThreshold;
if (!PyArg_ParseTuple(args, "O!f|O", &PyList_Type, &pyDictList, &userThreshold, &encodeBool)) {
    PyErr_SetString(PyExc_TypeError, "parameters must be a list of dictionnaries and a distance value.");
    return NULL;
}
bool bEncode;
if (encodeBool != NULL) {
if (PyObject_IsTrue(encodeBool))
    bEncode = true;
    Py_XDECREF(encodeBool);
}
Py_ssize_t nStructPairs = PyList_Size(pyDictList);
#ifdef DEBUG
    PySys_WriteStdout("Unpacking %d pairs of structure coordinates [contact distance is %f]\n", (int)nStructPairs, userThreshold);
#endif

    PyObject *PyListResults =  PyList_New(nStructPairs);
    PyObject *pStructAsDictRec, *pStructAsDictLig;

    atom_t *atomListRec, *atomListLig;
    int nAtomsRec, nAtomsLig;
    PyObject *PyTupleBuffer;

    for (int iStructPair = 0; iStructPair < (int)nStructPairs ; iStructPair++) {
#ifdef DEBUG
        PySys_WriteStdout("____Structure_____ %d\n", iStructPair);
#endif

        PyTupleBuffer = PyList_GetItem(pyDictList, iStructPair);
        pStructAsDictRec = PyTuple_GetItem(PyTupleBuffer, 0);
        pStructAsDictLig = PyTuple_GetItem(PyTupleBuffer, 1);

#ifdef DEBUG
        PySys_WriteStdout("REF COUNT pStructAsDictRec :: is %d\nREF COUNT pStructAsDictLig :: is %d\n", Py_REFCNT(pStructAsDictRec), Py_REFCNT(pStructAsDictLig));
#endif

        atomListRec = structDictToAtoms(pStructAsDictRec, &nAtomsRec);
        atomListLig = structDictToAtoms(pStructAsDictLig, &nAtomsLig);

#ifdef DEBUG
        PySys_WriteStdout("USED__REF COUNT pStructAsDictRec :: is %d\nUSED__REF COUNT pStructAsDictLig :: is %d\n", Py_REFCNT(pStructAsDictRec), Py_REFCNT(pStructAsDictLig));
#endif
    ccmapView_t *ccmapView;
    Py_BEGIN_ALLOW_THREADS
    bool bAtomic = false;
    ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                    ? &atomicContactMap
                    : &residueContactMap;
    ccmapView = computeMap(atomListRec, nAtomsRec, atomListLig, nAtomsLig, userThreshold, bEncode);
    Py_END_ALLOW_THREADS

    PyObject *cValue = ccmapViewToPyObject(ccmapView, bEncode);
    PyList_SetItem(PyListResults, iStructPair, cValue);
    destroyAtomList(atomListRec, nAtomsRec);
    destroyAtomList(atomListLig, nAtomsLig);
    destroyCcmapView(ccmapView);

#ifdef DEBUG
        PySys_WriteStdout("ATOM_DESTROYED__REF COUNT pStructAsDictRec :: is %d\nATOM_DESTROYED__REF COUNT pStructAsDictLig :: is %d\n", Py_REFCNT(pStructAsDictRec), Py_REFCNT(pStructAsDictLig));
        PySys_WriteStdout("Freeing json C pointer safely\n");
#endif

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

   return PyListResults;
}
/*
    Python API exposed as cmap
*/
static PyObject *ccmap_compute_flex(PyObject* self, PyObject* args, PyObject* kwargs) {
#ifdef DEBUG
PySys_WriteStderr("ccmap_compute_flex");
#endif

static char *kwlist[] = {"", "y", "dist", "atomic", "encode", NULL};

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
int iLen;
int jLen = 0;
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
PyObject *rValue = ccmapViewToPyObject(ccmapView, bEncode);
// Cleaning 
destroyAtomList(iAtomList, iLen);
if(coorDictJ != NULL)
    destroyAtomList(jAtomList, jLen);
destroyCcmapView(ccmapView);

return rValue;
}

static int myextension_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int myextension_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

/* MODULE DECLARATION */

static PyMethodDef ccmapMethods[] = {
    //...
    //{"compute",  ccmap_compute_base, METH_VARARGS,
     //"Compute a residue or atomic contact map."},
     {"cmap",   (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_flex, METH_VARARGS | METH_KEYWORDS,
     "Compute a residue or atomic contact map."},
    {"lzmap",   (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_zdock_pose_list, METH_VARARGS | METH_KEYWORDS,
     "YOUPI."},
    {
    "lmap",  ccmap_compute_list, METH_VARARGS,
    "Compute a protein-protein interface residue contact map."
     },
     {
    "zmap",  ccmap_compute_zdock_pose, METH_VARARGS,
    "Compute a contact map from a zdock-like encoded pose:\n"\
    "Parameters\n"\
    "\t2-tuple of dictorizedStructures\n"\
    "\tcontact distance float\n"\
    "\t3-uple of euler angles, aka docking pose-specific rotation\n"\
    "\t3-uple of translation of one translation vector components, aka docking-pose specific translation\n"\
    "\t3-uple of translation of one translation vector components, aka receptor translation offset\n"\
    "\t3-uple of translation of one translation vector components, aka ligand translation offset\n"\
    "\tEncoding flag, Boolean\n"\
    "Return\n\tA string if Encoding flag is false, a list of residue numbers otherwise\n"
    },
    //...
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "ccmap",
        NULL,
        sizeof(struct module_state),
        ccmapMethods,
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


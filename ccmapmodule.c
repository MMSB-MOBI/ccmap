#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
#include "mesh_io.h"
#include "transform_mesh.h"
#include "encode.h"
#include "ccmapmodule_utils.h"
#include "ccmapmodule_allocation.h"

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

/* ------------------    NON ZDOCK related API ------------------




*/

/*
    Python API exposed as cmap
*/
static PyObject *ccmap_compute(PyObject* self, PyObject* args, PyObject* kwargs) {
static char *kwlist[] = {"", "", "d", "atomic", "encode", NULL};

PyObject *atomicBool = NULL;
PyObject *encodeBool = NULL;
PyObject *coorDictI  = NULL;
PyObject *coorDictJ  = NULL;
atom_t *iAtomList    = NULL;
atom_t *jAtomList    = NULL;
float userThreshold  = 4.5;
bool bAtomic,bEncode;
int iLen             = 0;
int jLen             = 0;

if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
                                "O!|O!fOO", kwlist,\
                                &PyDict_Type, &coorDictI, \
                                &PyDict_Type, &coorDictJ, \
                                &userThreshold, \
                                &atomicBool,
                                &encodeBool )) {
    PyErr_SetString(PyExc_TypeError, "Wrong parameters");
    return NULL;
}

// Parsing map type and encoding scheme
setBooleanFromParsing(atomicBool, &bAtomic);
setBooleanFromParsing(encodeBool, &bEncode);

// Parsing coordinates

iAtomList = structDictToAtoms(coorDictI, &iLen);
if(coorDictJ != NULL)
    jAtomList = structDictToAtoms(coorDictJ, &jLen);
/*
PySys_WriteStderr("Running with %s coordinates sets dist=%f bAtom:%s bEncode:%s",\
                                    jAtomList == NULL ? "one" : "two",\
                                    userThreshold,\
                                    bAtomic ? "true" : "false",\
                                    bEncode ? "true" : "false");
*/
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

/*

    Python API lmap
    TODO
    Coordinates paramters => List of one or two elements
    Move the Pythread up ?

*/
static PyObject *ccmap_compute_list(PyObject* self, PyObject* args, PyObject* kwargs)
{

static char *kwlist[] = {"", "", "d", "atomic", "encode", NULL};

PyObject *pyDictArray_I   = NULL;
PyObject *pyDictArray_J   = NULL;
PyObject *encodeBool      = NULL;
PyObject *atomicBool      = NULL;

PyObject *PyListResults   = NULL;
PyObject *pStructAsDict_I = NULL;
PyObject *pStructAsDict_J = NULL;

Py_ssize_t structFrameLen;

atom_t **atomListList_I   = NULL;
atom_t **atomListList_J   = NULL;
int *nAtomsList_I         = NULL;
int *nAtomsList_J         = NULL;

float userThreshold = 4.5;
bool bEncode;
bool dual = false;
bool bAtomic = false;
ccmapView_t **ccmapViewList;


if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|OfOO", kwlist,\
                    &pyDictArray_I, \
                    &pyDictArray_J,\
                    &userThreshold, &atomicBool, &encodeBool)) {
    PyErr_SetString(PyExc_TypeError, "parameters must be a list of dictionnaries and a distance value.");
    return NULL;
}
#ifdef PYMEM_CHECK
fprintf(stderr, "ccmap_compute_list MEMORY ENTRY SUMMARY:\n");
fprintf(stderr, "pyDictList:%zd\n", pyDictArray_I->ob_refcnt);
if(pyDictArray_J != NULL)
    fprintf(stderr, "pyDictList:%zd\n", pyDictArray_J->ob_refcnt);
#endif

if (!PyArray_Check(pyDictArray_I)) {
    PyErr_SetString(PyExc_TypeError, "First coordinates set is not iterable");
    return NULL;
}

dual = pyDictArray_J != NULL;
if(dual) {
    if (!PyArray_Check(pyDictArray_J)) {
        PyErr_SetString(PyExc_TypeError, "Optional coordinates set is not iterable");
        return NULL;
    }  
    if( !PyArray_Equal(pyDictArray_I, pyDictArray_J) ) {
        PyErr_SetString(PyExc_TypeError, "Coordinates lists must have same sizes");
        return NULL;
    }
}

setBooleanFromParsing(encodeBool, &bEncode);
setBooleanFromParsing(atomicBool, &bAtomic);

structFrameLen = PyArray_Size(pyDictArray_I);

ccmap_compute_list_allocate(&ccmapViewList, \
                            &atomListList_I, &nAtomsList_I, \
                            &atomListList_J, &nAtomsList_J, \
                            structFrameLen, dual);

    // We off load from threads the loading of coordinates 
    for (int iStructPair = 0 ; iStructPair < (int)structFrameLen ; iStructPair++) {
        pStructAsDict_I                 = PyArray_GetItem(pyDictArray_I, iStructPair);
        Py_INCREF(pStructAsDict_I);
        atomListList_I[iStructPair]     = structDictToAtoms(pStructAsDict_I, &nAtomsList_I[iStructPair]);    
        Py_DECREF(pStructAsDict_I);

        if(dual) {
            pStructAsDict_J             = PyArray_GetItem(pyDictArray_J, iStructPair);
            Py_INCREF(pStructAsDict_J);
            atomListList_J[iStructPair] = structDictToAtoms(pStructAsDict_J, &nAtomsList_J[iStructPair]);    
            Py_DECREF(pStructAsDict_J);
        }
    }

/*
We wont be using function current scope python object through GIL recovery, we dont Py_INCREF them
*/

Py_BEGIN_ALLOW_THREADS
for (int i = 0; i < (int)structFrameLen ; i++) {
    atom_t *optAtomList_J = dual ? atomListList_J[i] : NULL;
    int opt_nAtom_J = dual ? nAtomsList_J[i] : -1;
    ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                    ? &atomicContactMap
                    : &residueContactMap;
    ccmapViewList[i] = computeMap(atomListList_I[i], nAtomsList_I[i], \
                                  optAtomList_J, opt_nAtom_J, \
                                  userThreshold, bEncode);
}
Py_END_ALLOW_THREADS

PyListResults   = PyList_New(structFrameLen);
for (int i = 0; i < (int)structFrameLen ; i++) {
    PyObject *cValue = ccmapViewToPyObject(ccmapViewList[i], bEncode);
    PyList_SetItem(PyListResults, i, cValue); // Reference cValue is stolen
}

ccmap_compute_list_cleanOnExit(ccmapViewList, \
                                atomListList_I, nAtomsList_I,\
                                atomListList_J, nAtomsList_J,\
                                structFrameLen, dual);

#ifdef PYMEM_CHECK
fprintf(stderr, "ccmap_compute_list MEMORY EXIT SUMMARY:\n");
fprintf(stderr, "pyDictArray_I:%zd\n", pyDictArray_I->ob_refcnt);
if(pyDictArray_J != NULL)
    fprintf(stderr, "pyDictArray_J:%zd\n", pyDictArray_J->ob_refcnt);

/*Garbage collactable on GIL release*/
if (pStructAsDict_I != NULL)
    fprintf(stderr, "pStructAsDict_I:%zd\n", pStructAsDict_I->ob_refcnt);
if (pStructAsDict_J != NULL)
    fprintf(stderr, "pStructAsDict_J:%zd (Borrowed pyDictList )\n", pStructAsDict_J->ob_refcnt);

fprintf(stderr, "PyListResults:%zd\n", PyListResults->ob_refcnt);

#endif
return PyListResults;
}

/*
Python API : "zmap"
*/

static PyObject *ccmap_compute_zdock_pose(PyObject *self, PyObject *args, PyObject* kwargs) 
{
static char *kwlist[] = {"", "", "", "", "offsetRec", "offsetLig", "d", "encode", "atomic", "apply", NULL};
PyObject *eulerArray         = NULL;
PyObject *translationArray   = NULL;
PyObject *offsetRecArray     = NULL;
PyObject *offsetLigArray     = NULL;
PyObject *encodeBool         = NULL;
PyObject *applyBool          = NULL;
PyObject *atomicBool         = NULL;

double *offsetLigVector      = NULL;
double *offsetRecVector      = NULL;
double *translation          = NULL;
double *eulerAngle           = NULL;
float userThreshold          = 4.5;

PyObject *pyDictRec          = NULL;
PyObject *pyDictLig          = NULL;

bool bEncode, bAtomic, bApply;
atom_t *atomListRec          = NULL; 
atom_t *atomListLig          = NULL;
int nAtomsRec, nAtomsLig;

ccmapView_t *ccmapView;

if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
                     "O!O!OO|OOfOOO", kwlist, \
                    &PyDict_Type, &pyDictRec,   \
                    &PyDict_Type, &pyDictLig,   \
                                  &eulerArray, \
                                  &translationArray,\
                                  &offsetRecArray,   \
                                  &offsetLigArray,   \
                    &userThreshold,
                    &encodeBool,                   
                    &atomicBool,
                    &applyBool               )) {
    PyErr_SetString(PyExc_TypeError, "Wrong parameters");
    return NULL;
}

setBooleanFromParsing(encodeBool, &bEncode);
setBooleanFromParsing(applyBool,  &bApply);
setBooleanFromParsing(atomicBool, &bAtomic);


atomListRec       = structDictToAtoms(pyDictRec, &nAtomsRec);
atomListLig       = structDictToAtoms(pyDictLig, &nAtomsLig);

eulerAngle        = unpackVector3(eulerArray); // euler angle to generate pose
translation       = unpackVector3(translationArray) ; // translation vector to generate pose
offsetLigVector   = unpackVector3(offsetLigArray); // translation to center ligand
offsetRecVector   = unpackVector3(offsetRecArray); // translation to center receptor

if (eulerAngle == NULL)
        PyErr_SetString(PyExc_TypeError, "Fail to unpack euler triplet");
if (translation == NULL) 
        PyErr_SetString(PyExc_TypeError, "Fail to unpack translation triplet");        
if (offsetLigVector == NULL) 
        PyErr_SetString(PyExc_TypeError, "Fail to unpack ligand offset triplet");
if (offsetRecVector == NULL)
        PyErr_SetString(PyExc_TypeError, "Fail to unpack receptor offset triplet");
if (offsetRecVector == NULL || offsetLigVector == NULL || translation == NULL || eulerAngle == NULL) {
    destroyVector3(eulerAngle);
    destroyVector3(translation);
    destroyVector3(offsetLigVector);
    destroyVector3(offsetRecVector);
}    

Py_BEGIN_ALLOW_THREADS
transformAtomList(atomListRec, NULL, offsetRecVector);

transformAtomList(atomListLig, eulerAngle, offsetLigVector);
transformAtomList(atomListLig, NULL, translation);

ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                ? &atomicContactMap
                : &residueContactMap;
ccmapView = computeMap(atomListRec, nAtomsRec, atomListLig, nAtomsLig, userThreshold, bEncode);
Py_END_ALLOW_THREADS

PyObject *rValue = ccmapViewToPyObject(ccmapView, bEncode);

if (bApply) {
    if ( ! backMapCoordinates(atomListRec, pyDictRec) ) {
        PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
        return NULL;
    }
    if ( ! backMapCoordinates(atomListLig, pyDictLig) ) {
        PyErr_SetString(PyExc_TypeError, "coordinates backmapping failed\n");
        return NULL;
    }
}

destroyVector3(eulerAngle);
destroyVector3(translation);
destroyVector3(offsetRecVector);
destroyVector3(offsetLigVector);
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
static char *kwlist[] = {"", "", "", "", "offsetRec", "offsetLig", "d", "encode", "atomic", NULL};
PyObject *pyDictRec       = NULL;
PyObject *pyDictLig       = NULL;
PyObject *eulerArrayArray = NULL;
PyObject *transArrayArray = NULL;
PyObject *offsetRecArray  = NULL;
PyObject *offsetLigArray  = NULL;
PyObject *encodeBool      = NULL;
PyObject *atomicBool      = NULL;
Py_ssize_t nPose;

double *offsetLigVector   = NULL;
double *offsetRecVector   = NULL;
double **eulerTriplets    = NULL;
double **transTriplets    = NULL;

ccmapView_t **ccmapViews;
bool bEncode, bAtomic;

float userThreshold       = 4.5;

int nAtomsRec             = 0;
int nAtomsLig             = 0;
atom_t *atomListRec       = NULL;
atom_t *atomListLig       = NULL;
atom_t *atomListLigBuffer = NULL;

if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
                     "O!O!OO|OOfOO", kwlist, \
                    &PyDict_Type, &pyDictRec,   \
                    &PyDict_Type, &pyDictLig,   \
                                  &eulerArrayArray,   \
                                  &transArrayArray,   \
                                  &offsetRecArray,   \
                                  &offsetLigArray,   \
                    &userThreshold,
                    &encodeBool,
                    &atomicBool                )) {
    PyErr_SetString(PyExc_TypeError, "Wrong parameters");
    return NULL;
}
setBooleanFromParsing(encodeBool, &bEncode);
setBooleanFromParsing(atomicBool, &bAtomic);

if( !PyArray_Equal(eulerArrayArray, transArrayArray) ) {
    PyErr_SetString(PyExc_TypeError, "Transformations tuples lists are of unequal length");
    return NULL;
}
if (offsetLigArray != NULL) {  
    offsetLigVector =  unpackVector3(offsetLigArray);    
    if (offsetLigVector == NULL ) {
        PySys_WriteStderr("Fail to unpack ligand offset triplet");
        return NULL;
    }
}
if (offsetRecArray != NULL) {  
    offsetRecVector =  unpackVector3(offsetRecArray);    
    if (offsetRecVector == NULL ) {
        destroyVector3(offsetLigVector);
        PySys_WriteStderr("Fail to unpack receptor offset triplet");
        return NULL;
    }
}


eulerTriplets = createListVector3(eulerArrayArray, &nPose);
if(eulerTriplets == NULL) {
    PySys_WriteStderr("Fail to unpack eulers tuples list");
    destroyVector3(offsetLigVector);
    destroyVector3(offsetRecVector);
    return NULL;
}
transTriplets = createListVector3(transArrayArray, &nPose);
if(transTriplets == NULL) {
    PySys_WriteStderr("Fail to unpack tranlsations tuples list");
    destroyVector3(offsetLigVector);
    destroyVector3(offsetRecVector);
    destroyListVector3(eulerTriplets, nPose);
    return NULL;
}

// Reading provided ligand, receptor conformations
atomListRec = structDictToAtoms(pyDictRec, &nAtomsRec);
atomListLig = structDictToAtoms(pyDictLig, &nAtomsLig);
atomListLigBuffer = structDictToAtoms(pyDictLig, &nAtomsLig);

// Alloc for results structures
ccmapViews = PyMem_New(ccmapView_t *, nPose); // MAybe not safe has access and subsequent malloc occur in thread below
ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool) = bAtomic\
                    ? &atomicContactMap
                    : &residueContactMap;

Py_BEGIN_ALLOW_THREADS
transformAtomList(atomListRec, NULL, offsetRecVector);
// Loop through all poses transformations
for (int iPose = 0 ; iPose < nPose ; iPose++) {
    // Transform Ligand
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
     {"cmap",   (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute, METH_VARARGS | METH_KEYWORDS,
        "Compute a residue or atomic contact map\n"\
        "One positional mandatory parameter:\n"\
        "\tOne structure as \"dictorized\" coordinates\n"\
        "One positional optional parameter:\n"\
        "\tAnother structure as \"dictorized\" coordinates\n"\
        "Three optional parameters:\n"
        "\t\"d\", contact distance float\n"\
        "\t\"encode\", Encoding flag, boolean. Default=False\n"\
        "\t\"atomic\", boolean. If True compute atomic contact map else compute residue contact map. Default=False\n"\
        "\nReturn\n\tA string if Encoding flag is false, a list of list of residue/atom numbers otherwise\n"
    },
    {
    "lcmap",  (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_list, METH_VARARGS | METH_KEYWORDS,
        "Compute a list of single or of pair or proteins residue or atomic contact maps\n"\
        "One positional mandatory parameter:\n"\
        "\tOne list of structures as \"dictorized\" coordinates\n"\
        "One positional optional parameter:\n"\
        "\tAnother list of structures as \"dictorized\" coordinates\n"\
        "Three optional parameters:\n"
        "\t\"d\", contact distance float\n"\
        "\t\"encode\", Encoding flag, boolean. Default=False\n"\
        "\t\"atomic\", boolean. If True compute atomic contact map else compute residue contact map. Default=False\n"\
    "\nReturn\n\tA string if Encoding flag is false, a list of list of residue/atom numbers otherwise\n"
    },
    {
    "zmap",  (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_zdock_pose, METH_VARARGS | METH_KEYWORDS,
    "Compute a contact map from a zdock-like encoded pose\n"\
    "WARNING: setting \"apply\" argument to True will modify the \"dictorized\" coordinates arguments!\n"\
    "Four positional mandatory parameters:\n"\
    "\tReceptor structure \"dictorized\" coordinates\n"\
    "\tLigand structure \"dictorized\" coordinates\n"\
    "\tList/Tuple of the euler angles, aka docking pose-specific rotation\n"\
    "\tList/Tuple of the translation vector, aka docking-pose specific translation\n"\
    "Five optional parameters:\n"
    "\t\"d\", contact distance float\n"\
    "\t\"offsetRec\", list/tuple of a translation vector, aka receptor translation offset\n"\
    "\t\"offsetLig\", list/tuple of a translation vector, aka ligand translation offset\n"\
    "\t\"encode\", Encoding flag, boolean. Default=False\n"\
    "\t\"atomic\", boolean. If True compute atomic contact map else compute residue contact map. Default=False\n"\
    "\t\"apply\", boolean. Apply the transformation to the provided receptor and ligand \"dictorized\" coordinates, default=False\n"\
    "\nReturn\n\tA string if Encoding flag is false, a list of residue/atom numbers otherwise\n"
    },
    {"lzmap",   (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_zdock_pose_list, METH_VARARGS | METH_KEYWORDS,
     "Compute the contact maps from a list of zdock-like encoded pose\n"\
    "Four positional mandatory parameters:\n"\
    "\tReceptor structure \"dictorized\" coordinates\n"\
    "\tLigand structure \"dictorized\" coordinates\n"\
    "\tList/Tuple of many list/3-uple of euler angles, aka docking pose-specific rotation\n"\
    "\tList/Tuple of many list/3-uple translation vectors, aka docking-pose specific translation\n"\
    "Four optional parameters:\n"
    "\t\"d\", contact distance float\n"\
    "\t\"offsetRec\", list/tuple of a translation vector, aka receptor translation offset\n"\
    "\t\"offsetLig\", list/tuple of a translation vector, aka ligand translation offset\n"\
    "\t\"encode\", Encoding flag, boolean. Default=False\n"\
    "\t\"atomic\", boolean. If True compute atomic contact map else compute residue contact map. Default=False\n"\
    "\nReturn\n\tA string if Encoding flag is false, a list of list of residue/atom numbers otherwise\n"
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


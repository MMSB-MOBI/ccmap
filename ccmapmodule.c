#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
#include "mesh_io.h"
#include "transform_mesh.h"
#include "encode.h"
#include "ccmapmodule_utils.h"

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
static PyObject *ccmap_compute_flex(PyObject* self, PyObject* args, PyObject* kwargs) {
#ifdef DEBUG
PySys_WriteStderr("ccmap_compute_flex");
#endif

static char *kwlist[] = {"", "y", "dist", "atomic", "encode", NULL};

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
jAtomList = NULL;
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

/*

    Python API lmap

*/
static PyObject *ccmap_compute_list(PyObject *self, PyObject *args)
{

PyObject *pyDictList       = NULL;
PyObject *encodeBool       = NULL;

PyObject *PyListResults    = NULL;
PyObject *pStructAsDictRec = NULL;
PyObject *pStructAsDictLig = NULL;
PyObject *PyTupleBuffer    = NULL;
Py_ssize_t nStructPairs;

atom_t *atomListRec        = NULL;
atom_t *atomListLig        = NULL;
int nAtomsRec, nAtomsLig;

float userThreshold = 4.5;
bool bEncode;

ccmapView_t *ccmapView;

if (!PyArg_ParseTuple(args, "O!f|O", &PyList_Type, &pyDictList, &userThreshold, &encodeBool)) {
    PyErr_SetString(PyExc_TypeError, "parameters must be a list of dictionnaries and a distance value.");
    return NULL;
}
setBooleanFromParsing(encodeBool, &bEncode);

nStructPairs = PyList_Size(pyDictList);
#ifdef DEBUG
    PySys_WriteStdout("Unpacking %d pairs of structure coordinates [contact distance is %f]\n", (int)nStructPairs, userThreshold);
#endif

PyListResults = PyList_New(nStructPairs);


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
Python API : "zmap"
*/

static PyObject *ccmap_compute_zdock_pose(PyObject *self, PyObject *args, PyObject* kwargs) 
{
static char *kwlist[] = {"", "", "", "", "offsetRec", "offsetLig", "distance", "encode", "apply", NULL};
PyObject *eulerTuple         = NULL;
PyObject *translationTuple   = NULL;
PyObject *offsetRec          = NULL;
PyObject *offsetLig          = NULL;
PyObject *encodeBool         = NULL;
PyObject *applyBool          = NULL;

double *offsetLigVector      = NULL;
double *offsetRecVector      = NULL;
double *translation          = NULL;
double *eulerAngle           = NULL;
float userThreshold          = 4.5;

PyObject *pyDictRec          = NULL;
PyObject *pyDictLig          = NULL;
bool bEncode, bApply;
bool bAtomic                 = false; // No atomic map in zdock context, yet
atom_t *atomListRec          = NULL; 
atom_t *atomListLig          = NULL;
int nAtomsRec, nAtomsLig;

ccmapView_t *ccmapView;

if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
                     "O!O!O!O!|O!O!fOO", kwlist, \
                    &PyDict_Type, &pyDictRec,   \
                    &PyDict_Type, &pyDictLig,   \
                    &PyTuple_Type, &eulerTuple, \
                    &PyTuple_Type, &translationTuple,\
                    &PyTuple_Type,  &offsetRec,   \
                    &PyTuple_Type,  &offsetLig,   \
                    &userThreshold,
                    &encodeBool,
                    &applyBool                )) {
    PyErr_SetString(PyExc_TypeError, "Wrong parameters");
    return NULL;
}

setBooleanFromParsing(encodeBool, &bEncode);
setBooleanFromParsing(applyBool, &bApply);

atomListRec       = structDictToAtoms(pyDictRec, &nAtomsRec);
atomListLig       = structDictToAtoms(pyDictLig, &nAtomsLig);

eulerAngle        = unpackVector3(eulerTuple); // euler angle to generate pose
translation       = unpackVector3(translationTuple) ; // translation vector to generate pose
offsetLigVector   = unpackVector3(offsetLig); // translation to center ligand
offsetRecVector   = unpackVector3(offsetRec); // translation to center receptor

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

#ifdef DEBUG
PySys_WriteStdout("Unpacking euler %.2f %.2f %.2f||translation vectors %.2f %.2f %.2f\n", \
    eulerAngle[0], eulerAngle[1], eulerAngle[2], \
    translation[0], translation[1], translation[2]);
//  PySys_WriteStdout("Unpacking %d rotation matrices [contact distance is %f]\n", (int)nRotMatrix, userThreshold);
#endif


Py_BEGIN_ALLOW_THREADS
/*
char titi[200];
stringifyAtom(&atomListRec[0], titi);
fprintf(stderr, "Z1:%s\n", titi);
*/
transformAtomList(atomListRec, NULL, offsetRecVector);
/*
stringifyAtom(&atomListRec[0], titi);
fprintf(stderr, "Z2:%s\n", titi);
*/
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
static char *kwlist[] = {"", "", "", "", "offsetRec", "offsetLig", "distance", "encode", NULL};
PyObject *pyDictRec       = NULL;
PyObject *pyDictLig       = NULL;
PyObject *eulerList       = NULL;
PyObject *transList       = NULL;
PyObject *offsetRec       = NULL;
PyObject *offsetLig       = NULL;
PyObject *encodeBool      = NULL;
Py_ssize_t nPose;

double *offsetLigVector   = NULL;
double *offsetRecVector   = NULL;
double **eulerTriplets    = NULL;
double **transTriplets    = NULL;

ccmapView_t **ccmapViews;
bool bEncode;
bool bAtomic              = false; // No atomic map in zdock context, yet
float userThreshold       = 4.5;

int nAtomsRec             = 0;
int nAtomsLig             = 0;
atom_t *atomListRec       = NULL;
atom_t *atomListLig       = NULL;
atom_t *atomListLigBuffer = NULL;

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
setBooleanFromParsing(encodeBool, &bEncode);

if( PyList_Size(eulerList) != PyList_Size(transList) ) {
    PyErr_SetString(PyExc_TypeError, "Transformations tuples lists are of unequal length");
    return NULL;
}
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


eulerTriplets = createListVector3(eulerList, &nPose);
if(eulerTriplets == NULL) {
    PySys_WriteStderr("Fail to unpack eulers tuples list");
    destroyVector3(offsetLigVector);
    destroyVector3(offsetRecVector);
    return NULL;
}
transTriplets = createListVector3(transList, &nPose);
if(transTriplets == NULL) {
    PySys_WriteStderr("Fail to unpack tranlsations tuples list");
    destroyVector3(offsetLigVector);
    destroyVector3(offsetRecVector);
    destroyListVector3(eulerTriplets, nPose);
    return NULL;
}

for (int k = 0; k < nPose ; k++) {
    PySys_WriteStdout("Unpacking %f %f %f // %f %f %f\n", \
                        eulerTriplets[k][0], eulerTriplets[k][1], eulerTriplets[k][2],\
                        transTriplets[k][0], transTriplets[k][1], transTriplets[k][2]\
                        );
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
     {"cmap",   (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_flex, METH_VARARGS | METH_KEYWORDS,
     "Compute a residue or atomic contact map. DO SOME DOC"},
    {
    "lmap",  ccmap_compute_list, METH_VARARGS,
        "Compute a list of protein-protein interface residue contact map."
     },
    {"lzmap",   (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_zdock_pose_list, METH_VARARGS | METH_KEYWORDS,
     "Compute the contact maps from a list of zdock-like encoded pose:\n"\
    "Four positional parameters:\n"\
    "\tReceptor structure \"dictorized\" coordinates\n"\
    "\tLigand structure \"dictorized\" coordinates\n"\
    "\tA list of 3-uple of euler angles, aka docking pose-specific rotation\n"\
    "\tA list of 3-uple of translation of one translation vector components, aka docking-pose specific translation\n"\
    "Four optional parameters:\n"
    "\t\"distance\", contact distance float\n"\
    "\t\"offsetRec\", tuple representation of translation vector components, aka receptor translation offset\n"\
    "\t\"offsetLig\", tuple representation of translation vector components, aka ligand translation offset\n"\
    "\t\"encode\", Encoding flag, boolean\n"\
    "\nReturn\n\tA string if Encoding flag is false, a list of list of residue numbers otherwise\n"\
    "\nExample:\n"
    },
     {
    "zmap",  (PyCFunction/*PyCFunctionWithKeywords*/)ccmap_compute_zdock_pose, METH_VARARGS | METH_KEYWORDS,
    "Compute a contact map from a zdock-like encoded pose:\n"\
    "Four positional parameters:\n"\
    "\tReceptor structure \"dictorized\" coordinates\n"\
    "\tLigand structure \"dictorized\" coordinates\n"\
    "\t3-uple of euler angles, aka docking pose-specific rotation\n"\
    "\t3-uple of translation of one translation vector components, aka docking-pose specific translation\n"\
    "Four optional parameters:\n"
    "\t\"distance\", contact distance float\n"\
    "\t\"offsetRec\", tuple representation of translation vector components, aka receptor translation offset\n"\
    "\t\"offsetLig\", tuple representation of translation vector components, aka ligand translation offset\n"\
    "\t\"encode\", Encoding flag, boolean\n"\
    "\t\"apply\", boolean. Apply the transformation to the provided receptor and ligand \"dictorized\" coordinates, default=False\n"\
    "\nReturn\n\tA string if Encoding flag is false, a list of residue numbers otherwise\n"\
    "\nExample:\n"
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


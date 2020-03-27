#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
// --------------------- Utility Functions ---------------------  XREF SANITY ?

// Returns a representaion of ccmapVie as a PyList if bEncode, PyString otherwise
PyObject *ccmapViewToPyObject(ccmapView_t *ccmapView, bool bEncode);
PyObject *ccmapViewsToPyObject(ccmapView_t **, int, bool);
int PyObject_AsDouble(PyObject *py_obj, double *x);
int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size);
int backMapCoordinates(atom_t *atomListRoot,  PyObject *pyDictObject);

double **createListVector3(PyObject *pyObject_List, Py_ssize_t *len);
double *unpackVector3(PyObject *pyObject_Tuple);
double **destroyListVector3(double **vList, Py_ssize_t len);
double *destroyVector3(double *vector);

int unpackChainID(PyObject *pListChainID, char **buffer);
int unpackString(PyObject *pListOfStrings, char ***buffer);
int unpackCoordinates(PyObject *pListCoor, double **buffer);
void freeBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n);
atom_t *structDictToAtoms(PyObject *pyDictObject, int *nAtoms);

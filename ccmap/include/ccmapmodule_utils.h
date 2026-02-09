#ifndef CCMAP_UTIL_H_INCLUDED
#define CCMAP_UTIL_H_INCLUDED
#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
//#include "numpy_headers.h"

// --------------------- Utility Functions ---------------------  XREF SANITY ?
bool PyArray_Equal(PyObject *arrayI, PyObject *arrayJ) ;
// Returns a representaion of ccmapVie as a PyList if bEncode, PyString otherwise
PyObject *ccmapViewToPyObject(ccmapView_t *ccmapView, bool bEncode);
PyObject *ccmapViewsToPyObject(ccmapView_t **, int, bool);
PyObject *ccmapViewToSasaPyDict(ccmapView_t *ccmapView);
void cmapViewAppendToSasaFrame(ccmapView_t *ccmapView, sasaFrame_t *sasaFrame, int iFrame);

int PyObject_AsDouble(PyObject *py_obj, double *x);
int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size);
int backMapCoordinates(atom_t *atomListRoot,  PyObject *pyDictObject);

PyObject *MyPyArray_GetItem(PyObject *pyObject_array, Py_ssize_t position);
bool MyPyArray_Check(PyObject *pyObject_arrayMaybe);
Py_ssize_t MyPyArray_Size(PyObject *pyObject_array);
double **createListVector3(PyObject *pyObject_List, Py_ssize_t *len);
double *unpackVector3(PyObject *pyObject_Tuple);
double **destroyListVector3(double **vList, Py_ssize_t len);
double *destroyVector3(double *vector);

int unpackChainID(PyObject *pListChainID, char **buffer);
int unpackString(PyObject *pListOfStrings, char ***buffer);
double *unpackCoordinates(PyObject *pListCoor);
void freeBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n);
atom_t *structDictToAtoms(PyObject *pyDictObject, int *nAtoms, float probeRadius, atom_map_t *aMap, int sasaHiRes);
void setBooleanFromParsing(PyObject *, bool *);

PyObject *buildPyValueSasaFrame(sasaFrame_t *sasaFrame);
atom_map_t *dictRadiiToAtomMapper(PyObject *atomRadiiPyDict);
//void PyObject_ToChar(PyObject *source, char *target);

#endif
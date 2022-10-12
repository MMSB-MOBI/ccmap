//File: python_utils.h
#ifndef PYTHON_UTILS_H
#define PYTHON_UTILS_H

#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include "miscellaneous.h"


bool PyArray_Equal(PyObject *arrayI, PyObject *arrayJ);
int PyObject_AsDouble(PyObject *py_obj, double *x);
int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size);
PyObject *MyPyArray_GetItem(PyObject *pyObject_array, Py_ssize_t position);
bool MyPyArray_Check(PyObject *pyObject_arrayMaybe);

Py_ssize_t MyPyArray_Size(PyObject *pyObject_array);

double **createListVector3(PyObject *pyObject_array, Py_ssize_t *len);
double **destroyListVector3(double **vList, Py_ssize_t len);
double *unpackVector3(PyObject *pyObject);
double *destroyVector3(double *vector);
int unpackChainID(PyObject *pListChainID, char **buffer);
int unpackString(PyObject *pListOfStrings, char ***buffer);
double *unpackCoordinates(PyObject *pListCoor);
void PyObject_ToChar(PyObject *source, char *target);
#endif

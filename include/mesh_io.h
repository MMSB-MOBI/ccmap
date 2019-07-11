//int unpackCoordinates(PyObject *pListCoor, double *buffer);
int PyObject_AsDouble(PyObject *py_obj, double *x);
int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size);

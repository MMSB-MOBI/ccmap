#include "python_utils.h"

bool PyArray_Equal(PyObject *arrayI, PyObject *arrayJ) {
    Py_ssize_t (*PyArray_SizeI)(PyObject *);
    Py_ssize_t (*PyArray_SizeJ)(PyObject *);

    PyArray_SizeI = PyList_Check(arrayI) ? &PyList_Size : &PyTuple_Size;
    PyArray_SizeJ = PyList_Check(arrayJ) ? &PyList_Size : &PyTuple_Size;
  
    return PyArray_SizeI == PyArray_SizeJ;
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

/* Get item at position in a tuple or a list */
PyObject *MyPyArray_GetItem(PyObject *pyObject_array, Py_ssize_t position) {
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

bool MyPyArray_Check(PyObject *pyObject_arrayMaybe) {
    if (PyList_Check(pyObject_arrayMaybe) )
        return true;
    if (PyTuple_Check(pyObject_arrayMaybe) )
        return true;
    
    return false;
}

Py_ssize_t MyPyArray_Size(PyObject *pyObject_array) {
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
    const char *s;
    Py_ssize_t sLen;

    for (i = 0 ; i < n ; i++) {
        pItem = PyList_GetItem(pListChainID, i);
        Py_INCREF(pItem);
        if(!PyUnicode_Check(pItem))
            return NULL;
        
        s = PyUnicode_AsUTF8AndSize(pItem, &sLen);
        Py_DECREF(pItem);
        
        (*buffer)[i] = s[0];       
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

    const char* s;
    Py_ssize_t sLen;
    for (i = 0; i < n ; i++) {
        pItem = PyList_GetItem(pListOfStrings, i);
        Py_INCREF(pItem);

        if(!PyUnicode_Check(pItem))
            return NULL;
        
        s = PyUnicode_AsUTF8AndSize(pItem, &sLen);
        Py_DECREF(pItem);
     
        (*buffer)[i] = PyMem_New(char, sLen + 1);
        for (int j = 0 ; j < sLen; j++) 
            (*buffer)[i][j] = s[j];
        
        (*buffer)[i][sLen] = '\0';
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



/* FKN Unicode new style api 
    PyListObject *Py_keylist = PyDict_Keys(pyDictObject);
    Py_ssize_t n_key         = PyList_Size(Py_keylist);
    PyObject *buf;
    const char *toto;
    Py_ssize_t sLen;
    for(int i = 0 ; i < n_key ; i++) {
        buf = PyList_GetItem(Py_keylist, i);
        if(PyUnicode_Check(buf)){
            fprintf(stderr, "key %d is unicode of length %d\n", i,(int)PyUnicode_GET_LENGTH(buf));
            toto = PyUnicode_AsUTF8AndSize(buf, &sLen);
            fprintf(stderr, "key is %s\n", toto);
        }
    }
*/
void PyObject_ToChar(PyObject *source, char *target) {
    PyObject *pyStr = PyUnicode_AsUTF8String(source);
    const char *s = PyBytes_AS_STRING(pyStr);
    strcpy(target, s); 
    
    Py_DECREF(pyStr);
}
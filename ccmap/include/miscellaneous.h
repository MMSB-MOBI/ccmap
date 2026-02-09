//File: miscellaneous.h
#ifndef MISCELLANEOUS_H
#define MISCELLANEOUS_H

#include "math.h"
#include "stdio.h"

/*
#ifdef AS_PYTHON_EXTENSION
#include <Python.h>
#endif
*/

#if DEBUG
    void printOnContextStderr(char *debugCharArray);
#endif

double euclideanDistance3(float x0, float y0, float z0, float x1, float y1, float z1);

void strip(char *dest, char *src);

typedef enum { false, true } bool;

#endif

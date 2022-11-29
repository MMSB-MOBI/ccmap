//File: mesh_default.h
#ifndef MESH_DEFAULT_H
#define MESH_DEFAULT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_string.h"
#include "miscellaneous.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

//Utility functions
int concatenate(char **dest, char *src);
int popChar(char **dest, char c);

#endif

//File: parameters.h
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

void stringToThreeFloats(char *input, float (*vector)[3]);
void parseTransform(char *eulerString, char *translateString, float (*eulers)[3], float (*translate)[3]);

#endif

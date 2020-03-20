//File: default.h
#ifndef DEFAULT_H
#define DEFAULT_H

typedef enum { false, true } bool;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

//Utility functions
int concatenate(char **dest, char *src);
int popChar(char **dest, char c);

#endif
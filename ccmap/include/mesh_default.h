//File: mesh_default.h
#ifndef MESH_DEFAULT_H
#define MESH_DEFAULT_H

typedef enum { false, true } bool;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>

typedef struct string string_t;

#define STRING_CHUNCK (1024)//(2)

struct string {
    int length;
    size_t capacity;
    char *value;
    bool (*append)(string_t *, char *);
    void (*dump)(string_t *);
    string_t *(*copy)(string_t *);
    char*(*toChar)(string_t *);
};

//Utility functions
int concatenate(char **dest, char *src);
int popChar(char **dest, char c);

bool appendString(string_t *this, char *word);
string_t *copyString(string_t *this);
string_t *createString(void); 
string_t *destroyString(string_t *);
void dumpString(string_t *);
char *toCharString(string_t *);

#endif

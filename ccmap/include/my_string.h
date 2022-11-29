//File: mesh_default.h
#ifndef MY_STRING_H
#define MY_STRING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include "miscellaneous.h"

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

typedef struct stringList {
    string_t **elem;
    int len;
} stringList_t;



bool appendString(string_t *this, char *word);
string_t *copyString(string_t *this);
string_t *createString(void); 
string_t *destroyString(string_t *);
void dumpString(string_t *);
char *toCharString(string_t *);

stringList_t *splitAndCreateStringList(char *input, char sep);
stringList_t *destroyStringList(stringList_t *);
bool compareStringToChar(string_t *self, char *other, bool bStrip);
#endif

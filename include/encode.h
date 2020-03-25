//File: encode.h
#ifndef ENCODE_H
#define ENCODE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "molecular_object.h"

int chainLen(residue_t *ResidueList);
int contactIndex(int index1, int index2, int max2);
int *encodeContactMap(residue_t *ResidueList, int lenLigList, int lenRecList, unsigned int *finalLen);

void printTable(int *ContactList, unsigned int len);
int *copyTable(int *table, int lenTable );
#endif

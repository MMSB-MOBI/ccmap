//File: encode.h
#ifndef ENCODE_H
#define ENCODE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define TABLE_CHUNCK_SZ (1000)
#define ENCODE_IJ2K(iRec, jLig, nLig) ( (iRec) * (nLig) + (jLig) )
#include "molecular_object.h"

int chainLen(residue_t *ResidueList);
int contactIndex(int index1, int index2, int max2);
//int *encodeContactMap(residue_t *ResidueList, int lenLigList, int lenRecList, unsigned int *finalLen);
unsigned int *encodeContactMap(residue_t *iResidueList, residue_t *jResidueList, unsigned int *totalContacts);
void printTable(int *ContactList, unsigned int len);
unsigned int *copyTable(unsigned int *table, int lenTable );
#endif

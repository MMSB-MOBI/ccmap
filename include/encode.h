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

unsigned int *encodeContactMapResidue(residue_t *iResidueList, residue_t *jResidueList, unsigned int *totalContacts);
unsigned int *encodeContactMapAtomic(atom_t *iatomList, atom_t *jAtomList, atomPair_t *ccList, unsigned int *totalContacts);
unsigned int *copyTable(unsigned int *table, int lenTable );

#endif

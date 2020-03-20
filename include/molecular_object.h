//File: molecular_object.h
#ifndef MOLECULAR_OBJECT_H
#define MOLECULAR_OBJECT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pdb_coordinates.h"

typedef struct residue {
    struct atom *elements;
    int nAtoms;
    int index;
    char *resID;
    char chainID;
    char *resName;
    struct residue *nextResidueList;
    struct residue *prevResidueList;
    struct residue **contactResidueList; // list of residue ptr dynmically incremented
    int nContacts;
   // struct residue
} residue_t;

typedef struct atom {
    struct residue *belongsTo;
    struct cell *inCell;
    char chainID;
    char *resID;
    float x;
    float y;
    float z;
    char *resName;
    char *name;
    struct atom *nextAtomList;
    struct atom *nextCellAtom;
    struct atom *nextResidueAtom;
} atom_t;

typedef struct atomPair {
    struct atom *a;
    struct atom *b;
    double dist;
    struct atomPair *next;
} atomPair_t;


// Display molecular objects content
void printResidueList(residue_t *residueList);
void printResidue(residue_t *residue);
void printAtomList(atom_t *atomList);
void stringifyAtom(atom_t *atom, char *atomString);
void stringifyResidue(residue_t *residue, char *residueString);
void jsonifyResidue(residue_t *residue, char *jsonString);
void jsonArrayifyAtom(atom_t *atom, char *atomString, bool bCoordinates);

residue_t *createResidue(atom_t *atom, int n);
residue_t *createResidueList(atom_t * atomList);
void fuseResidueLists(residue_t *iResidueList, residue_t *jResidueList);
residue_t *destroyResidueList(residue_t *residueList);
residue_t *destroyResidue(residue_t *residue);
atom_t *destroyAtomList(atom_t *atomList, int nAtom);
atom_t *destroyAtom(atom_t *atom);

int CreateAtomListFromPdbContainer(pdbCoordinateContainer_t *pdbCoordinateContainer, atom_t *atomList);
atom_t *readFromArrays(int nAtoms, double *x, double *y, double *z, char *chainID, char **resID, char **resName, char **name);
void freeAtomListCreatorBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n);

atom_t *legacy_readCoordinates(char *fname, int *_nAtom);
#endif
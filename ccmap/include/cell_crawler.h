//File: cell_crawler.h
#ifndef CELL_CRAWLER_H
#define CELL_CRAWLER_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "molecular_object.h"
#include "mesh_default.h"
#include "cell.h"

#ifdef AS_PYTHON_EXTENSION
#include <Python.h>
#endif

#define ATOM_CC_LIST_CHUNCK (1000)


typedef struct updaterStruct updaterStruct_t;
typedef struct cellCrawler cellCrawler_t;

// Defines the type of operation to perform when we find
// pairwise atomic contact
struct updaterStruct {
    uint64_t maxSize;
    atomPair_t *atomContactList;
    uint64_t totalByAtom;
    uint64_t totalByResidue;
    bool (*updaterFn)(cellCrawler_t*, atom_t*, atom_t*, double);
};

// Defines the kind of atom pairwise distance to compute
// At a given cell
// Hook the updater to call for dist below threshold
struct cellCrawler {
    void (*enumerator)(cellCrawler_t *, cell_t *, cell_t *);
    bool dual;
    updaterStruct_t *updater;
    double threshold;
    bool (*atomPairProcess)(cellCrawler_t*, atom_t*, atom_t*);
};

bool processPairwiseDistance(cellCrawler_t*, atom_t*, atom_t*);

cellCrawler_t *createCellCrawler(bool atomic, bool dual, double dist);
void extendCellCrawler(cellCrawler_t *cellCrawler);
cellCrawler_t *destroyCellCrawler(cellCrawler_t *cellCrawler);

void pairwiseCellEnumerate(cellCrawler_t *, cell_t *refCell, cell_t *targetCell);
void pairwiseCellEnumerateDual(cellCrawler_t *, cell_t *refCell, cell_t *targetCell);

// Atomic lvl contact registring
bool updateAtomContact(cellCrawler_t *, atom_t *iAtom, atom_t *jAtom, double dist);
// Residue lvl contact registring
bool updateResidueContact(cellCrawler_t *, atom_t *iAtom, atom_t *jAtom, double dist);

double distance(atom_t *iAtom, atom_t *jAtom);

#endif

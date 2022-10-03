//File: sasa.h
#ifndef SASA_H
#define SASA_H

#include "molecular_object.h"
#include "chem_constants.h"
#include "fibonacci.h"

typedef struct residue_sasa {
    //residue_t *residue; bound to residue lifetime we dont want that
    char resname[81];
    int res_index;   // The index(numbering/rank) of the residue
    char resID[81];
    char chainID;
    float buried;
    float nominal;
    float frac;
} residue_sasa_t;

typedef struct sasaResults {
    residue_sasa_t *residueSasaList;
    uint16_t length;
} sasaResults_t;

sasaResults_t *computeSasaResults(residueList_t *residueList);
sasaResults_t *destroySasaResults(sasaResults_t *sasaResults);

string_t      *jsonifySasaResults(sasaResults_t *sasaResults);
#endif

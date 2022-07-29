//File: sasa.h
#ifndef SASA_H
#define SASA_H

#include "molecular_object.h"
#include "chem_constants.h"
#include "fibonacci.h"

typedef struct residue_sasa {
    residue_t *residue;
    float buried;
    float nominal;
    float frac;
} residue_sasa_t;

typedef struct sasaResults {
    residue_sasa_t *residueSasaList;
    u_int16_t length;
} sasaResults_t;

sasaResults_t *computeSasaResults(residueList_t *residueList);
sasaResults_t *destroySasaResults(sasaResults_t *sasaResults);

string_t      *jsonifySasaResults(sasaResults_t *sasaResults);


#endif

//File: sasa.h
#ifndef SASA_H
#define SASA_H

#include "molecular_object.h"
#include "chem_constants.h"

typedef struct residue_sasa {
    residue_t *residue;
    float sasa;
} residue_sasa_t;

typedef struct sasaResults {
    residue_sasa_t *residueSasaList;

} sasaResults_t;

sasaResults_t *sasaCore(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, float probeRadius);

typedef struct cellSasaCrawler {
    // TO DO
} cellSasaCrawler_t;

#endif

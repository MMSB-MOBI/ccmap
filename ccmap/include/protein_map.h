//File: protein_map.h
#ifndef PROTEIN_MAP_H
#define PROTEIN_MAP_H

#include "molecular_object.h"

typedef struct atom_info {
        int atype;              // TYPES OF ATOMS
        int nb;                 // NUMBER OF ATOMS IN RESIDUES
        float vrad;             // VAN DER WAALS RADII
        int keyresidue;         // RESIDUE KEY FROM ENUM PROTEIN MAP
        int keyatom;            // ATOM KEY FROM ENUM PROTEIN MAP
} atom_info_t;

typedef float (*VDW_MAPPER)(atom_t*);

float getVdwRadius(atom_t*);

#endif

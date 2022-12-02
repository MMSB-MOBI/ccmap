//File: convex.h
#ifndef CONVEX_H
#define CONVEX_H

#include "molecular_object.h"
#include "mesh.h"

void buildSurfaces(meshContainer_t *meshContainer);
int buildSphere(atom_t *atom, cell_t *cell, meshContainer_t *meshContainer);
bool surfaceExplorerPredicate(cell_t *cell);
#endif

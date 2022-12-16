//File: convex.h
#ifndef CONVEX_H
#define CONVEX_H

#include "molecular_object.h"
#include "mesh.h"

bool buildSurfaces(meshContainer_t *meshContainer, bool force);
int buildSphere(atom_t *atom, cell_t *cell, meshContainer_t *meshContainer);
bool surfaceExplorerPredicate(cell_t *cell);
int voxelEvaluate(cell_t *currCell, cell_t *centerCell, double norm);
int voxelEvaluateCartesian(meshContainer_t *meshContainer, cell_t *currCell, atom_t *centerAtom);
int countCorners(meshContainer_t *meshContainer, atom_t *atom, cell_t *cell);
setCells_t *getSurfaceCells(atom_t *atom, meshContainer_t *meshContainer);
#endif

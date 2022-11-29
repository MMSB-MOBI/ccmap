//File: pathfinder.h
#ifndef PATHFINDER_H
#define PATHFINDER_H
#include "molecular_object.h"
#include <math.h>
#include <mesh.h>

#define DEFAULT_BWFS 999999

typedef struct path {
    cell_t **cells;
    int len;
    cell_t *start;
    cell_t *stop;
} path_t;

typedef struct offsets {
    int i;
    int j;
    int k;
    float abs_dist;
} offsets_t;

path_t *searchForPath(meshContainer_t *meshContainer, atom_t *start, atom_t *stop);
short int sortNeighboursByMeshDistance(cell_t *currentCell, cell_t *goal, mesh_t *mesh, offsets_t moves[]);
int cmpOffsetfunc (const void * a, const void * b);
float c_dist(cell_t *a, cell_t *b);
void exploreCell(meshContainer_t *meshContainer, cell_t *currentCell, int nStepFromStart, cell_t *start_cell, cell_t *end_cell);
bool areSameCells(cell_t *a, cell_t *b);
path_t *backtrack(meshContainer_t *meshContainer, cell_t *startCell, cell_t *stopCell);
cell_t *walkBack(cell_t *currCell, cell_t *targetCell, meshContainer_t *meshContainer);
path_t *destroyPath(path_t *path);

#endif

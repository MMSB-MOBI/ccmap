//File: pathfinder.h
#ifndef PATHFINDER_H
#define PATHFINDER_H
#include "molecular_object.h"
#include "convex.h"
#include <math.h>
#include <mesh.h>

#define DEFAULT_BWFS 999999

typedef struct path {
    cell_t **cells;
    int len;
    cell_t *start;
    cell_t *stop;
    double patch_len_up;
    double patch_len_lo;
} path_t;

typedef struct offsets {
    int i;
    int j;
    int k;
    float abs_dist;
} offsets_t;

path_t *searchForPath(meshContainer_t *meshContainer,  char *type,\
                      atom_t *start, atom_t *stop);
short int sortNeighboursByMeshDistance(cell_t *currentCell, cell_t *goal, mesh_t *mesh, offsets_t moves[]);
int cmpOffsetfunc (const void * a, const void * b);
float c_dist(cell_t *a, cell_t *b);
void exploreCell(meshContainer_t *meshContainer, bool (*cellPredicate)(cell_t*),\
                 cell_t *currentCell, int nStepFromStart, cell_t *start_cell, cell_t *end_cell);
bool pointExplorerPredicate(cell_t *cell);
bool areSameCells(cell_t *a, cell_t *b);
path_t *backtrack(meshContainer_t *meshContainer, cell_t *startCell, cell_t *stopCell, char *type);
cell_t *walkBack(cell_t *currCell, cell_t *targetCell, meshContainer_t *meshContainer, bool (*cellPredicate)(cell_t*));
void reallocErrorLog(int a, int b, char type[]);

int createRecordArraysFromPath(path_t *self, meshContainer_t *meshContainer,\
                                double **x, double **y, double **z, char **chainID,\
                                char ***resID, char ***resName, char ***pearl_name, char segID, double spacing);    
       

path_t *destroyPath(path_t *path);

#endif

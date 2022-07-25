//File: fibonacci.h
#ifndef FIBONACCI_H
#define FIBONACCI_H

#include "mesh_default.h"

#define FIBO_K14 (377)
#define FIBO_K13 (233)


#include "chem_constants.h"

typedef struct spot {
    float x;
    float y;
    float z;
    bool buried;
} spot_t;

typedef struct fibo_grid {
    struct spot center;
    struct spot *spots;
    int n_spots;
    /* data */
} fibo_grid_t;

fibo_grid_t *computeFiboGrid(float x, float y, float z, float radius);

fibo_grid_t *createFiboGrid(int n);
fibo_grid_t *destroyFiboGrid(fibo_grid_t *fibo_grid);

void printFiboGrid(fibo_grid_t *fibo_grid);
string_t *jsonifyFiboGrid(fibo_grid_t *fibo_grid);

#endif

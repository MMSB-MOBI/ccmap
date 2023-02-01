//File: fibonacci.h
#ifndef FIBONACCI_H
#define FIBONACCI_H

#include "mesh_default.h"
#include "miscellaneous.h"

#define FIBO_K11 (89)
#define FIBO_K10 (55)

#define FIBO_K14 (377)
#define FIBO_K13 (233)

#define FIBO_K19 (4181)
#define FIBO_K20 (6765)


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
    float radius;
    /* data */
} fibo_grid_t;

fibo_grid_t *computeFiboGrid(float x, float y, float z, float radius, int resolutionLevel);

fibo_grid_t *createFiboGrid(int n);
fibo_grid_t *destroyFiboGrid(fibo_grid_t *fibo_grid);

void printFiboGrid(fibo_grid_t *fibo_grid);
string_t *jsonifyFiboGrid(fibo_grid_t *fibo_grid);
void FiboSpherePairProcess(fibo_grid_t *iFiboGrid, fibo_grid_t *jFiboGrid);
void computeFiboSphereASA(fibo_grid_t *iFiboGrid, float *totalSurface, float *buriedSurface);
void updateFiboGrid(fibo_grid_t *fibo_grid, float dx, float dy, float dz);
int generateGridPointCartesian(fibo_grid_t *self, double *x, double *y, double *z, bool surfOnly);
#endif

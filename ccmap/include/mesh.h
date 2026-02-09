//File: mesh.h
#ifndef MESH_H
#define MESH_H

#include "mesh_default.h"
#include "molecular_object.h"
#include "cell_crawler.h"
#include "cell.h"
#include "encode.h"
#include "sasa.h"

#ifdef AS_PYTHON_EXTENSION
#include <Python.h>
#endif

typedef struct ccmapView {
    char *asJSON;
    unsigned int *asENCODE;
    size_t encodeLen;
    sasaResults_t *sasaResults;
} ccmapView_t;

typedef struct ccmapResults {
    residueList_t *iResidueList;
    residueList_t *jResidueList;
    cellCrawler_t *cellCrawler;
    sasaResults_t *sasaResults;
    bool fused;
} ccmapResults_t;

typedef struct meshContainer {
    struct mesh *mesh;
    struct cell **filledCells;
    int nFilled;
    float step; // Grid step
    double x_min;
    double y_min;
    double z_min;
    int voxel_shell_thickness; //
    int nVoxels;
} meshContainer_t;

typedef struct mesh {
    int iMax;
    int jMax;
    int kMax;
    int n;
    cell_t ***grid;
} mesh_t;

typedef struct setCells {
    cell_t **cells;
    int size;
} setCells_t;

// ** API computing Fn **
ccmapView_t *residueContactMap(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, bool bEncoded, bool bASA);
ccmapView_t *atomicContactMap(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, bool bEncoded, bool bASA);

// Experimental SASA API
void generateSasaCloudToPdb(pdbCoordinateContainer_t **cloudPdbContainer, atom_t *iAtomList, int iAtom, double cellDim);


//Results display
char *jsonifyContactList(residueList_t *residueList); // BUMP TO string_t* TO DO
string_t *jsonifyAtomPairList(ccmapResults_t *ccmapResults);

// Mesh engine Fn
ccmapResults_t *ccmapCore(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, bool bAtomic, bool bASA);

void meshCrawler(meshContainer_t *meshContainer,  cellCrawler_t *cellCrawler);
cell_t ** vectorizeMesh(mesh_t *mesh);
// Constructors/Destructors
meshContainer_t *createMeshContainer(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double step);
mesh_t *createMesh(int iDim, int jDim, int kDim);
meshContainer_t *destroyMeshContainer(meshContainer_t *container);
mesh_t *destroyMesh(mesh_t *i_mesh);
ccmapResults_t *createCcmapResults(cellCrawler_t *, residueList_t *, residueList_t *);
ccmapResults_t *destroyCcmapResults (ccmapResults_t *results);
ccmapView_t *createCcmapView(void);
ccmapView_t *destroyCcmapView(ccmapView_t *);
setCells_t *destroySetCells(setCells_t*);
// Utilities
void getBoundariesCartesian_DUAL(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, atom_t *minCoor, atom_t *maxCoor);
void getBoundariesCartesian(atom_t * atomList, int nAtom, atom_t *minCoor, atom_t *maxCoor);
void cartesianToMesh(atom_t *atom, int *i, int *j, int *k, float step, atom_t minCoor, int vxsthickness);
void meshToCartesian(meshContainer_t *meshContainer, int i, int j, int k, double *x, double *y, double *z);
// kinda public version of meshToCartesian
cell_t *getCellFromAtom(meshContainer_t *meshContainer, atom_t *atom);
//void mesh(atom_t * atomList, int nAtom, double step);
float c_dist(cell_t *a, cell_t *b);
int manh_dist(cell_t *a, cell_t *b);
cell_t *selectFromSetCellByProx(setCells_t *set, cell_t *target, int mode);

//Debug Fn
void printContactList(residue_t *residueList);

void printMesh(mesh_t *mesh);
void dumpMeshContent(meshContainer_t *meshContainer);
void dumpCellContent(cell_t *cell);
void meshDummy(int a, int b, int c);
void printResidueCellProjection(char *resID, char chainID, meshContainer_t *meshContainer, residue_t *residueList);
void atomListInContact(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, int iAtomStatus[], int jAtomStatus[]);
void inspect(meshContainer_t *meshContainer, int _i, int _j, int _k, atom_t *iAtom);

void appendVoxelToPdbContainer(pdbCoordinateContainer_t *pdbContainer, meshContainer_t *meshContainer, char vID, bool buriedHidden);
#endif

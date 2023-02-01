//File: molecular_object.h
#ifndef MOLECULAR_OBJECT_H
#define MOLECULAR_OBJECT_H

#define  _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
//#include "numpy_headers.h"


#include "pdb_coordinates.h"
#include "mesh_default.h"
#include "fibonacci.h"
#include "chem_constants.h"
#include "atom_mapper.h"
#ifdef AS_PYTHON_EXTENSION
    #include "python_utils.h"


    #define NO_IMPORT_ARRAY
    #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
    #define PY_ARRAY_UNIQUE_SYMBOL CCMAP_ARRAY_API
    #include <Python.h>
    #include "numpy/arrayobject.h"
#endif

typedef struct residue {
    struct atom *elements;
    int nAtoms;
    int index;   // The index(numbering/rank) of the residue
    char *resID;
    char chainID;   
    char *ext_chainID; // MD analysis segid mutli char storage
    char *resName;
    struct residue *nextResidueList;
    struct residue *prevResidueList;
    struct residue **contactResidueList; // list of residue ptr dynmically incremented
    int nContacts;
   // struct residue
} residue_t;

typedef struct residueList {
    struct residue *root;
    uint16_t length;
} residueList_t;

typedef struct atom {
    struct residue *belongsTo;
    struct cell *inCell;
    char chainID;
    char *ext_chainID; // MD analysis segid mutli char storage
    char *resID;
    float x;
    float y;
    float z;
    char *resName;
    char *name;
    float _radiusASA;
    struct atom *nextAtomList;
    struct atom *nextCellAtom;
    struct atom *nextResidueAtom;
    unsigned int index;
    struct fibo_grid *f_grid;
    float _VDWradius; // for login purpose
} atom_t;

typedef struct atomPair {
    struct atom *a;
    struct atom *b;
    double dist;
    struct atomPair *next;
} atomPair_t;

typedef struct coordinates {
    float x;
    float y;
    float z;
   // atom_t *atom;
} coordinates_t;

typedef struct coorFrame {
    struct coordinates **coordinates;
    int nbFrame;
    int nbAtom;
} coorFrame_t;

typedef struct sasa2uple {
    double SASA;
    double frac;
} sasa2uple_t;

typedef struct sasaFrame {
    struct sasa2uple **sasa2upleArray;
    int nbFrame;
    int nbRes;
    residueList_t *residueList;
} sasaFrame_t;


// Display molecular objects content
void printResidueList(FILE *stream, residueList_t *residueList);
//unsigned int residueListLen(residue_t *);
void printResidue(FILE *stream, residue_t *residue);
void printAtomList(atom_t *atomList, FILE *stream);
void stringifyAtom(atom_t *atom, char *atomString);
void stringifyResidue(residue_t *residue, char *residueString);
void jsonifyResidue(residue_t *residue, char *jsonString);
void jsonArrayifyAtom(atom_t *atom, char *atomString, bool bCoordinates);
void jsonifyAtomPair(atomPair_t *atomPair, char *jsonString);

bool applyCoordinates(atom_t *atomListFrom, atom_t *atomListTo);
residue_t *createResidue(atom_t *atom, int n);
residueList_t *createResidueList(atom_t * atomList);
residueList_t *fuseResidueLists(residueList_t *iResidueList, residueList_t *jResidueList);
residueList_t *destroyResidueList(residueList_t *residueList);
residue_t *destroyResidue(residue_t *residue);
atom_t *destroyAtomList(atom_t *atomList, int nAtom);
atom_t *destroyAtom(atom_t *atom);
void clearAtomFields(atom_t *atom);


unsigned int atomPairListLen(atomPair_t *);
unsigned int atomListLen(atom_t *);
// Following one should v its prtotype expanded and spreaded in code
atom_t *createBareboneAtom(int n, double x, double y, double z, char chainID, char *resID, \
                           char *resName, char *name);
void updateCoordinates(atom_t *atomList, coordinates_t *coordinates);
atom_t *CreateAtomListFromPdbContainer(pdbCoordinateContainer_t *pdbCoordinateContainer, int *nAtom, atom_map_t *aMap, float probeRadius, int sasaResLvl);
atom_t *readFromArrays(int nAtoms, double *x, double *y, double *z, char *chainID, char **resID, char **resName, char **name, atom_map_t *aMap, float probeRadius, int sasaHiRes);
void freeAtomListCreatorBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n);

sasaFrame_t *createSasaFrame(atom_t *atomList, int nbFrame);
sasaFrame_t *destroySasaFrame(sasaFrame_t *sasaFrame);
coorFrame_t *destroyCoorFrame(coorFrame_t *coorFrame, int optIndex);

int appendFiboGridToPdbContainer(pdbCoordinateContainer_t *cloudPdbContainer,\
                                atom_t *atomList, int nbAtom, char segID, bool hiddenBuried);

atom_t *legacy_readCoordinates(char *fname, int *_nAtom);
#ifdef AS_PYTHON_EXTENSION
    atom_t *readFromNumpyArrays(PyArrayObject *positions, PyArrayObject *names,\
                                PyArrayObject *resnames,  PyArrayObject *resids,\
                                PyArrayObject *segids, atom_map_t *aMap, float probeRadius, int resolutionLevel);
    coorFrame_t *createFrameFromPyArrayList(PyObject *positionArrayList);
    atom_t *readFromNumpyArraysFrame(coorFrame_t **coorFrame, PyObject *positionFrame,\
                                PyArrayObject *_names, PyArrayObject *_resnames,\
                                PyArrayObject *_resids, PyArrayObject *_segids,\
                                atom_map_t *aMap, float probeRadius, int resolutionLevel);

#endif

#endif


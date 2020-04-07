#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "decoygen.h"
#include "mesh.h"

/*
typedef struct atom {
    struct residue *belongsTo;
    struct cell *inCell;
    char chainID;
    char *resID;
    float x;
    float y;
    float z;
    char *resName;
    char *name;
    struct atom *nextAtomList;
    struct atom *nextCellAtom;
    struct atom *nextResidueAtom;
} atom_t;
*/



/*
    A COPY OF pdb_ccordinates.c: transformPdbCoordinateContainer,
    adapted to iterate over a list of mesh atom_t pointers
*/
void transformAtomList(atom_t *atomListRoot, double euler[3], double translation[3]) {
    atom_t *atomCurrent = atomListRoot;

    float nX, nY, nZ;

    if (translation != NULL) {
    /*    #ifdef AS_PYTHON_EXTENSION
        #ifdef DEBUG
            PySys_WriteStdout("Applying translation %g, %g, %g\n", translation[0], translation[1], translation[2]);
        #endif
        #endif
    */
        while (atomCurrent != NULL) {
    /*
           #ifdef DEBUG
            float oX = atomCurrent->x;
            float oY = atomCurrent->y;
            float oZ = atomCurrent->z;
            #endif
    */
            atomCurrent->x += translation[0];
            atomCurrent->y += translation[1];
            atomCurrent->z += translation[2];
        /*
        #ifdef DEBUG
            #ifdef AS_PYTHON_EXTENSION
                PySys_WriteStdout("[trans]X: %f -> %f, Y:%f -> %f, Z : %f->%f\n", oX, atomCurrent->x, oY, atomCurrent->y, oZ, atomCurrent->z);
            #else
            fprintf(stderr,"[trans]X: %f -> %f, Y:%f -> %f, Z : %f->%f\n", oX, atomCurrent->x, oY, atomCurrent->y, oZ, atomCurrent->z);
            #endif
            #endif
        */
            atomCurrent = atomCurrent->nextAtomList;
        }
    }
    atomCurrent = atomListRoot;
    if (euler != NULL) {
    /*
        #ifdef AS_PYTHON_EXTENSION
        #ifdef DEBUG
            PySys_WriteStdout("Applying euler rotations %g, %g, %g\n", euler[0], euler[1], euler[2]);
        #endif
        #endif
    */
        while (atomCurrent != NULL) {
            rotateAtom(atomCurrent->x,  atomCurrent->y,  atomCurrent->z,
                       &nX    ,  &nY    ,  &nZ    ,
                      euler[0],  euler[1], euler[2]);
            
    /*
            #ifdef AS_PYTHON_EXTENSION
            #ifdef DEBUG
                PySys_WriteStdout("[euler]X: %f -> %f, Y:%f -> %f, Z : %f->%f\n", atomCurrent->x, nX, atomCurrent->y, nY, atomCurrent->z, nZ);
            #endif
            #endif
    */
            atomCurrent->x = nX;
            atomCurrent->y = nY;
            atomCurrent->z = nZ;
            atomCurrent = atomCurrent->nextAtomList;
        }
    }
}


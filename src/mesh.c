#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mesh.h"
#ifndef ENCODE
#include "encode.h"
#define ENCODE
#endif
#ifdef AS_PYTHON_EXTENSION
#include <Python.h>
#endif




// Pop the last element of a string if it matches provided char c
// The length of the string is effectively reduced by one, this is the returned value
int popChar(char **dest, char c) {
    int i = 0;
    while( (*dest)[i] != '\0' ) {
        i++;
    }
    if ( (*dest)[i - 1] == c ) {
        (*dest)[i - 1] = '\0';
         *dest = realloc( *dest, (i) * sizeof(char) );
    }
    return i;
}
/*
strcpy Copies the C string pointed by source into the array
pointed by destination, including the terminating null character
(and stopping at that point).
strlen(const char *str) computes the length
of the string str up to, but not including the terminating null character.
*/



// Both chains are supposed to have a teminating null char Alternatively dest can be a NULL string pointer
// Concatenate two strings in the 1st one, return the length of the string not including the '\0' termination character
int concatenate(char **dest, char *src) {
    int oldStringLength = *dest != NULL ? strlen(*dest) : 0;
    // oldStrinLength is the index at which we will start to append/copy the src srting into dest string
    int new_buf_size = oldStringLength + strlen(src) + 1; // Both string length + slot for '\0'

    *dest = realloc( *dest, new_buf_size * sizeof(char) );
    strcpy(&((*dest)[oldStringLength]), src);
    if ((*dest)[new_buf_size - 1] != '\0') {
        printf("String copy buffer termination error");
    }
    // This line should not be necessary
    (*dest)[new_buf_size - 1] = '\0';

    return new_buf_size - 1;
}



void dumpBuffer(char *buffer, int bufSize) {
    printf("checking buffer\n");
    for (int i = 0; i < bufSize ; i++) {
        printf("[%d]'%c'", i, buffer[i]);
        if (buffer[i] == '\0') {
            printf("<----HERE\n");
        } else {
            printf("\n");
        }
    }
}
// To communicate w/ python process
// Return a JSON with residue index and contact
// [ [1,2], [4,6] .. ]

char *residueContactMap(atom_t * atomList, int nAtom, double ctc_dist) {
    residue_t *residueList = createResidueList(atomList);
#ifdef AS_PYTHON_EXTENSION
PySys_WriteStdout("Using residueContactMap function \n");
#endif

#ifdef DEBUG
#ifdef AS_PYTHON_EXTENSION
        PySys_WriteStdout("Computing residue contact map w/ %.2g Angstrom step\n", ctc_dist);
#endif
    printf("Computing residue contact map w/ %.2g Angstrom step\n", ctc_dist);
#endif

    double step = ctc_dist;
    meshContainer_t *results = createMeshContainer(atomList, nAtom, NULL, 0, step);
    int nPairs;
    enumerate(results, ctc_dist, &nPairs, false);

    char *jsonString = jsonifyContactList(residueList);

#ifdef DEBUG
    printContactList(residueList);
    printf("%s\n", jsonString);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStderr("%s\n", jsonString);
#endif
#endif
    residueList = destroyResidueList(residueList);
    results = destroyMeshContainer(results);
    return jsonString;
}
// Debugging function to list the cell coordinates of a specified residues projected atoms
static void printResidueCellProjection(char *resID, char chainID, meshContainer_t *meshContainer, residue_t *residueList) {
    char atomString[81];

    residue_t *residuePtr = residueList;
    atom_t *atomPtr = NULL;
    while(residuePtr != NULL) {
        if( (strcmp(resID, residuePtr->resID) == 0) && residuePtr->chainID == chainID ) {
            printf("CellProj:: \"%s\" \"%c\"\n", residuePtr->resID, residuePtr->chainID);
            atomPtr = residuePtr->elements;
            while(atomPtr != NULL) {
                stringifyAtom(atomPtr, atomString);
                printf("CellProj:: %s [%d, %d, %d]\n", atomString, atomPtr->inCell->i, atomPtr->inCell->j, atomPtr->inCell->k);
                atomPtr = atomPtr->nextResidueAtom;
            }
        }
        residuePtr = residuePtr->nextResidueList;
    }
}

void atomListInContact(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double step, int iAtomStatus[], int jAtomStatus[]) {

    meshContainer_t *results = createMeshContainer(iAtomList, iAtom, jAtomList, jAtom, step);

    cell_t *cur_cell;
    for (int i = 0; i < iAtom ; i++) {
        iAtomStatus[i] = 0;
        cur_cell = iAtomList[i].inCell;
        if (cur_cell->jMemberCount > 0)
            iAtomStatus[i] = 1;
    }

    for (int j = 0; j < jAtom ; j++) {
        jAtomStatus[j] = 0;
        cur_cell = jAtomList[j].inCell;
        if (cur_cell->iMemberCount > 0)
            jAtomStatus[j] = 1;
    }

    results = destroyMeshContainer(results);
}

int *residueContactMap_DUAL(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, unsigned int *finalLen) {
    residue_t *iResidueList = createResidueList(iAtomList);
    residue_t *jResidueList = createResidueList(jAtomList);

#ifdef DEBUG
#ifdef AS_PYTHON_EXTENSION
        PySys_WriteStdout("Computing residue contact map w/ %.2g Angstrom step\n", ctc_dist);
#endif
    printf("Computing residue contact map w/ %.2g Angstrom step\n", ctc_dist);
#endif
    double step = ctc_dist;
    meshContainer_t *results = createMeshContainer(iAtomList, iAtom, jAtomList, jAtom, step);
    /* Inspecting atom projection */
    // 101_B_CE1 and 121_1_OE1 cell coordinates ?
    // printResidueCellProjection(" 101", 'B', results, iResidueList);
    //printResidueCellProjection(" 121", 'A', results, jResidueList);
    int nPairs;
    enumerate(results, ctc_dist, &nPairs, true);

    // Link the two residues list
    int jlen=chainLen(jResidueList);
    int ilen=chainLen(iResidueList);
    fuseResidueLists(iResidueList, jResidueList);
    int *ccmap=encodeContactMap(iResidueList, jlen, ilen, finalLen);
    /*printf("Number of contacts : %d\n", *finalLen);*/

    // We dont jsonify anymore, relying on int encoding by Julia
    //char *jsonString = jsonifyContactList(iResidueList);

    #ifdef DEBUG
        printContactList(iResidueList);
        printf("%s\n", jsonString);
    #ifdef AS_PYTHON_EXTENSION
        PySys_WriteStderr("%s\n", jsonString);
        PySys_WriteStderr("\n" );
    #endif
    #endif


    iResidueList = destroyResidueList(iResidueList);
    results = destroyMeshContainer(results);
    return ccmap;
}

// connect j with trailing i element
// Both are chained list , fuse them should be safe in terms of destruction
void fuseResidueLists(residue_t *iResidueList, residue_t *jResidueList) {
    residue_t *head = iResidueList;
    while(head->nextResidueList != NULL)
        head = head->nextResidueList;
    head->nextResidueList = jResidueList;
    jResidueList->prevResidueList = head;
}

void printContactList(residue_t *residueList) {
    //residueList = iterate over residue list, iterate over its contact -> stringify residue pair;
    residue_t *residue_curr = residueList;
    residue_t *residue_partner = NULL;
    char residueString_current[81];
    char residueString_partner[81];

    while(residue_curr != NULL) {
        stringifyResidue(residue_curr, residueString_current);
        for (int i = 0 ; i <  residue_curr->nContacts ; i++) {
            residue_partner = residue_curr->contactResidueList[i];
            stringifyResidue(residue_partner, residueString_partner);
#ifdef DEBUG
            printf("%s <> %s\n", residueString_current, residueString_partner);
#endif
        }
        residue_curr = residue_curr->nextResidueList;
    }
}

char *jsonifyContactList(residue_t *residueList) {
    //residueList = iterate over residue list, iterate over its contact jsonIfy;
    residue_t *residue_curr = residueList;
    residue_t *residue_partner = NULL;
    char residueJsonString_current[81];
    char residueJsonString_partner[81];
    char jsonHeader[] = "{\"type\":\"contactList\", \"data\":[\0";

    char jsonRootElemTag[] = "{\"root\":\0";
    char jsonPartnerElemTagOpen[] = ",\"partners\":[\0";
    char jsonPartnerElemTagClose[] = "]}\0";
    char jsonFooter[] = "]}\0";

    int bufSize = (strlen(jsonHeader) + 1);
    char *jsonStringTotal = malloc( bufSize * sizeof(char) );
    strcpy(jsonStringTotal, jsonHeader);

#ifdef DEBUG
    printf("Starting jsonification\n");
#endif
    //int lastCharPosition = strlen(jsonStringTotal) - 1;
    while(residue_curr != NULL) {
        if(residue_curr->nContacts > 0) {
            concatenate(&jsonStringTotal, jsonRootElemTag);
            jsonifyResidue(residue_curr, residueJsonString_current);
            concatenate(&jsonStringTotal, residueJsonString_current);
            concatenate(&jsonStringTotal, jsonPartnerElemTagOpen);

            for (int i = 0 ; i <  residue_curr->nContacts ; i++) {
                residue_partner = residue_curr->contactResidueList[i];
                jsonifyResidue(residue_partner, residueJsonString_partner);
                concatenate(&jsonStringTotal, residueJsonString_partner);
                if (i < residue_curr->nContacts - 1)
                    concatenate(&jsonStringTotal, ",\0");
            }
            concatenate(&jsonStringTotal, jsonPartnerElemTagClose);
            if (residue_curr->nextResidueList != NULL)
                concatenate(&jsonStringTotal, ",\0");
        }
        residue_curr = residue_curr->nextResidueList;
    }
    popChar(&jsonStringTotal, ',');
    concatenate(&jsonStringTotal, jsonFooter);
#ifdef DEBUG
    printf("jsonification performed successfully\n");
#endif
    return jsonStringTotal;
}

void jsonifyResidue(residue_t *residue, char *jsonString) {
    sprintf(jsonString, "{\"resID\" : \"%s\", \"chainID\":\"%c\"}", residue->resID, residue->chainID);
}


residue_t *createResidue(atom_t *atom, int n) {
    residue_t *residue = malloc(sizeof(residue_t));
    residue->nAtoms = n;
    residue->nContacts = 0;
    residue->resName = malloc( (strlen(atom->resName) + 1) * sizeof(char) );
    strcpy(residue->resName, atom->resName);
    residue->resID = malloc( (strlen(atom->resID) + 1) * sizeof(char) );
    strcpy(residue->resID, atom->resID);
    residue->chainID = atom->chainID;
    residue->elements = atom;
    residue->elements->nextResidueAtom = NULL;
    residue->contactResidueList = NULL;
    residue->index = n;
    residue->prevResidueList = NULL;
    residue->nextResidueList = NULL;

    return residue;
}


residue_t *createResidueList(atom_t * atomList) {
    residue_t *residue_root = NULL;
    residue_t *residue_head = NULL;
    atom_t *curr_atom = atomList;//&atomList[0];
    atom_t *prev_atom = NULL;
    int nResidue = 0;

    residue_root = createResidue(curr_atom, nResidue);
    residue_head = residue_root;

    curr_atom->belongsTo = residue_head;
    prev_atom = curr_atom;
    curr_atom = curr_atom->nextAtomList;

    while(curr_atom != NULL) {
        if (strcmp(curr_atom->resID, residue_head->resID) == 0
            && curr_atom->chainID == residue_head->chainID) {
            prev_atom->nextResidueAtom = curr_atom;
            residue_head->nAtoms += 1;
        } else {
            nResidue++;
            residue_head->nextResidueList = createResidue(curr_atom, nResidue);
            residue_head->nextResidueList->prevResidueList = residue_head;
            residue_head = residue_head->nextResidueList;
        }
        curr_atom->nextResidueAtom = NULL;
        prev_atom = curr_atom;
        curr_atom->belongsTo = residue_head;

        curr_atom = curr_atom->nextAtomList;
    }
#ifdef DEBUG
    printf("Created %d residues\n", nResidue);
#endif
    return residue_root;
}

// MEMORY ALLOCATION OF atom LIST
atom_t *readFromArrays(int nAtoms, double *x, double *y, double *z, char *chainID, char **resID, char **resName, char **name) {
    atom_t *atomList = malloc(nAtoms * sizeof(atom_t));
    for (int n = 0 ; n < nAtoms ; n++) {
        atomList[n].nextAtomList = NULL;
        atomList[n].belongsTo = NULL;
        atomList[n].inCell = NULL;
        atomList[n].nextResidueAtom = NULL;

        atomList[n].resID = malloc( (strlen(resID[n]) + 1) * sizeof(char));
        strcpy(atomList[n].resID, resID[n]);

        atomList[n].resName = malloc( (strlen(resName[n]) + 1) * sizeof(char));
        strcpy(atomList[n].resName, resName[n]);

        atomList[n].name = malloc( (strlen(name[n]) + 1) * sizeof(char));
        strcpy(atomList[n].name, name[n]);

        atomList[n].x = x[n];
        atomList[n].y = y[n];
        atomList[n].z = z[n];
        atomList[n].chainID = chainID[n];

        if (n > 0)
            atomList[n - 1].nextAtomList = &atomList[n];
    }
    return atomList;
}

void printResidueList(residue_t *residueList) {
    residue_t *curr_residue = residueList;
    printf("PRINTING RESIDUE LIST\n");

    printResidue(curr_residue);
    while(curr_residue->nextResidueList != NULL) {
        curr_residue = curr_residue->nextResidueList;
        printResidue(curr_residue);
    }
}



void stringifyAtom(atom_t *atom, char *atomString) {
    sprintf(atomString, "%s %s %s %c %g %g %g", atom->name, atom->resName, atom->resID, atom->chainID, atom->x, atom->y, atom->z);
}
void stringifyResidue(residue_t *residue, char *residueString) {
    sprintf(residueString, "%s %s %c [i=%d, n=%d]", residue->resName, residue->resID, residue->chainID, residue->index, residue->nAtoms);
}

void printAtomList(atom_t *atomList) {
    atom_t *curr_atom = &atomList[0];
    char atomString[81];
    while(curr_atom->nextAtomList != NULL) {
        stringifyAtom(curr_atom, atomString);
        printf("%s\n", atomString);
        curr_atom = curr_atom->nextAtomList;
    }
    stringifyAtom(curr_atom, atomString);
    printf("%s\n", atomString);
}

void printResidue(residue_t *residue) {
    atom_t *curr_atom;
    char atomString[81];
    char residueString[81];
    stringifyResidue(residue, residueString);
    printf("RES: %s (%d)\n", residueString, residue->nAtoms);
    curr_atom = residue->elements;
    while(curr_atom->nextResidueAtom != NULL) {
        stringifyAtom(curr_atom, atomString);
        printf("%s\n", atomString);
        curr_atom = curr_atom->nextAtomList;
    }
    stringifyAtom(curr_atom, atomString);
    printf("%s\n", atomString);
    printf("--------------------\n");
}


mesh_t *createMesh(int iDim, int jDim, int kDim) {
    mesh_t *i_mesh = malloc(sizeof(mesh_t));
    i_mesh->iMax = iDim;
    i_mesh->jMax = jDim;
    i_mesh->kMax = kDim;
    i_mesh->n = iDim * jDim * kDim;
#ifdef DEBUG
    printf ("Creating a %d * %d * %d GRID\n", i_mesh->iMax, i_mesh->jMax, i_mesh->kMax);
#endif
    int n = 0;
    i_mesh->grid = malloc(i_mesh->iMax * sizeof(cell_t**));
    for ( int i = 0 ; i < i_mesh->iMax ; i++ ) {
        i_mesh->grid[i] = malloc(i_mesh->jMax * sizeof(cell_t*));
        for ( int j = 0 ; j < i_mesh->jMax ; j++ ) {
            i_mesh->grid[i][j] = malloc(i_mesh->kMax * sizeof(cell_t));
            for (int k = 0 ;  k < i_mesh->kMax ; k++) {
                i_mesh->grid[i][j][k].n = n++;
                i_mesh->grid[i][j][k].i = i;
                i_mesh->grid[i][j][k].j = j;
                i_mesh->grid[i][j][k].k = k;
                i_mesh->grid[i][j][k].memberCount = 0;
                i_mesh->grid[i][j][k].neighbourCount = 0;
                i_mesh->grid[i][j][k].members = NULL;

                i_mesh->grid[i][j][k].iMembers = NULL;
                i_mesh->grid[i][j][k].jMembers = NULL;

                i_mesh->grid[i][j][k].iMemberCount = 0;
                i_mesh->grid[i][j][k].jMemberCount = 0;
            }
        }
    }
#ifdef DEBUG
    printf("Mesh created\n");
#endif
    return i_mesh;
}

void printMesh(mesh_t *mesh) {
    printf("Mesh %dx%dx%d elements\n", mesh->iMax, mesh->jMax, mesh->kMax);
    for (int i = 0; i < mesh->iMax; i++) {
        for (int j = 0; j < mesh->jMax; j++) {
            for (int k = 0; k < mesh->kMax; k++) {
                printf ("%5d ", mesh->grid[i][j][k].n);
            }
            printf("\n");
        }
        printf("-------------------\n");
    }
}

// Go through none empty cells
// Get its following cells neighbours
void enumerate(meshContainer_t *meshContainer, double ctc_dist, int *nPairs, bool dualBool) {
    mesh_t *mesh = meshContainer->mesh;
    cell_t **cellList = meshContainer->filledCells;
    int nCells = meshContainer->nFilled;
    //int (*functionPtr)(int,int);
    void (*pCellEnumerator)(cell_t*, cell_t*, double , int*, int*);

    if (dualBool)
        pCellEnumerator = &pairwiseCellEnumerate_DUAL;
    else
        pCellEnumerator = &pairwiseCellEnumerate;

    cell_t ***grid = mesh->grid;
    cell_t *cur_cell;
    int kStart, jStart;
    bool extractBool = ctc_dist > 0.0 ? true : false;

    int nDist = 0;
    int nContacts = 0;
#ifdef DEBUG
    printf("Enumerating Distance between %d grid cells (Grid step is %g)\n", nCells, meshContainer->step);
#endif
    for (int c = 0 ; c < nCells ; c++) {
        cur_cell = cellList[c];
#ifdef DEBUG
        printf("Neighbours of cell %d (%d %d %d)\n", cur_cell->n, cur_cell->i, cur_cell->j, cur_cell->k);
#ifdef AS_PYTHON_EXTENSION
     //   PySys_WriteStderr("Neighbours of cell %d (%d %d %d)\n", cur_cell->n, cur_cell->i, cur_cell->j, cur_cell->k);
#endif
#endif
        // List Neighbouring cells and enumerate the pairwise atomic distances
        for (int i = cur_cell->i ; i <= cur_cell->i + 1 ; i++) {
            if (i >= mesh->iMax) break;
            jStart = i == cur_cell->i ? cur_cell->j : cur_cell->j - 1;
            for (int j = jStart ; j <= cur_cell->j + 1 ; j++) {
                if (j < 0) continue;
                if (j >= mesh->jMax) break;
                kStart = cur_cell->k - 1;
                if (i == cur_cell->i && j == cur_cell->j) kStart = cur_cell->k;
                for (int k = kStart ; k <= cur_cell->k + 1 ; k++) {
                    if (k < 0) continue;
                    if (k >= mesh->kMax) break;
                    if(!extractBool) {
                        printf("%d ", grid[i][j][k].n);
                        continue;
                    }
                    //parwiseCellEnumerate(cur_cell, &grid[i][j][k], ctc_dist, &nContacts, &nDist); // Enumerate contact w/ cell
                    (*pCellEnumerator)(cur_cell, &grid[i][j][k], ctc_dist, &nContacts, &nDist);
                }
            }
        }
    }
#ifdef DEBUG
    printf("\n ---> %d distances computed for a total of %d residue contacts\n", nDist, nContacts);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("\n ---> %d distances computed for a total of %d residue contacts\n", nDist, nContacts);
#endif
#endif
}

void pairwiseCellEnumerate(cell_t *refCell, cell_t *targetCell, double ctc_dist, int *nContacts, int *nDist) {
    atom_t *iAtom, *jAtom;

    char iAtomString[81];
    char jAtomString[81];

    if(refCell->memberCount == 0 || targetCell->memberCount == 0) return;
#ifdef DEBUG
    printf("\n*********\nPairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("Pairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#endif
#endif
    iAtom = refCell->members;
    while(iAtom != NULL) {

    #ifdef DEBUG
        stringifyAtom(iAtom, iAtomString);
    #endif
        jAtom = targetCell->members;
        while(jAtom != NULL) {
            if (jAtom != iAtom) {

            #ifdef DEBUG
                stringifyAtom(jAtom, jAtomString);
            #endif
                if(distance(iAtom, jAtom) < ctc_dist) {
                    (*nContacts) += updateContactList(iAtom, jAtom);
                }
                (*nDist) = (*nDist)+ 1;
#ifdef DEBUG
                printf("[[Dnum %d]] %s [%d %d %d]// %s [%d %d %d] ==> %.2g\n", *nDist, iAtomString, refCell->i, refCell->j, refCell->k, jAtomString, targetCell->i, targetCell->j, targetCell->k, distance(iAtom, jAtom));
#ifdef AS_PYTHON_EXTENSION
                PySys_WriteStdout("[[Dnum %d]] %s [%d %d %d]// %s [%d %d %d] ==> %.2g\n", *nDist, iAtomString, refCell->i, refCell->j, refCell->k, jAtomString, targetCell->i, targetCell->j, targetCell->k, distance(iAtom, jAtom));
#endif
#endif
            }
            jAtom = jAtom->nextCellAtom;
        }
        iAtom = iAtom->nextCellAtom;
    }
}

void pairwiseCellEnumerate_DUAL(cell_t *refCell, cell_t *targetCell, double ctc_dist, int *nContacts, int *nDist) {
    atom_t *iAtom, *jAtom;

    char iAtomString[81];
    char jAtomString[81];

    if(refCell->memberCount == 0 || targetCell->memberCount == 0) return;
#ifdef DEBUG
    printf("\n*********\nDUAL Pairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("DUAL Pairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#endif
#endif


    iAtom = refCell->iMembers;
    while(iAtom != NULL) {

    #ifdef DEBUG
        stringifyAtom(iAtom, iAtomString);
    #endif
        jAtom = targetCell->jMembers;
        while(jAtom != NULL) {
            if (jAtom != iAtom) {
            #ifdef DEBUG
                stringifyAtom(jAtom, jAtomString);
            #endif
                if(distance(iAtom, jAtom) < ctc_dist) {
                    (*nContacts) += updateContactList_DUAL(iAtom, jAtom);
                }
                (*nDist) = (*nDist)+ 1;
#ifdef DEBUG
                printf("DUAL [[Dnum %d]] %s [%d %d %d]// %s [%d %d %d] ==> %.2g\n", *nDist, iAtomString, refCell->i, refCell->j, refCell->k, jAtomString, targetCell->i, targetCell->j, targetCell->k, distance(iAtom, jAtom));
#ifdef AS_PYTHON_EXTENSION
                PySys_WriteStdout("DUAL [[Dnum %d]] %s [%d %d %d]// %s [%d %d %d] ==> %.2g\n", *nDist, iAtomString, refCell->i, refCell->j, refCell->k, jAtomString, targetCell->i, targetCell->j, targetCell->k, distance(iAtom, jAtom));
#endif
#endif
            }
            jAtom = jAtom->nextCellAtom;
        }
        iAtom = iAtom->nextCellAtom;
    }
// Reverse i/j members lookup
    iAtom = refCell->jMembers;
    while(iAtom != NULL) {

    #ifdef DEBUG
        stringifyAtom(iAtom, iAtomString);
    #endif
        jAtom = targetCell->iMembers;
        while(jAtom != NULL) {
            if (jAtom != iAtom) {
            #ifdef DEBUG
                stringifyAtom(jAtom, jAtomString);
            #endif
                if(distance(iAtom, jAtom) < ctc_dist) {
		  //(*nContacts) += updateContactList_DUAL(iAtom, jAtom);
		  // Reverse again for the contactList, to always have root residues from the receptor
		    (*nContacts) += updateContactList_DUAL(jAtom, iAtom);
                }
                (*nDist) = (*nDist)+ 1;
#ifdef DEBUG
                printf("DUAL [[Dnum %d]] %s [%d %d %d]// %s [%d %d %d] ==> %.2g\n", *nDist, iAtomString, refCell->i, refCell->j, refCell->k, jAtomString, targetCell->i, targetCell->j, targetCell->k, distance(iAtom, jAtom));
#ifdef AS_PYTHON_EXTENSION
                PySys_WriteStdout("DUAL [[Dnum %d]] %s [%d %d %d]// %s [%d %d %d] ==> %.2g\n", *nDist, iAtomString, refCell->i, refCell->j, refCell->k, jAtomString, targetCell->i, targetCell->j, targetCell->k, distance(iAtom, jAtom));
#endif
#endif
            }
            jAtom = jAtom->nextCellAtom;
        }
        iAtom = iAtom->nextCellAtom;
    }

}

int updateContactList(atom_t *iAtom, atom_t *jAtom){
    if(iAtom->belongsTo == jAtom->belongsTo) return 0;

    residue_t *iResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? iAtom->belongsTo : jAtom->belongsTo;
    residue_t *jResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? jAtom->belongsTo : iAtom->belongsTo;

    for (int i = 0; i < iResidue->nContacts; i++) {
        if(iResidue->contactResidueList[i] == jResidue) {
          //  printf("Contact already knwow between residues indexed %d,%d\n", iResidue->index, jResidue->index);
            return 0;
        }
    }
    //printf("ADDING a new contact between residues indexed %d,%d\n", iResidue->index, jResidue->index);
    iResidue->nContacts++;
    iResidue->contactResidueList = realloc( iResidue->contactResidueList, iResidue->nContacts * sizeof(residue_t*) );
    iResidue->contactResidueList[iResidue->nContacts - 1] = jResidue;
    return 1;
}


int updateContactList_DUAL(atom_t *iAtom, atom_t *jAtom){
    if(iAtom->belongsTo == jAtom->belongsTo) return 0;

    residue_t *iResidue = iAtom->belongsTo;
    residue_t *jResidue = jAtom->belongsTo;

    for (int i = 0; i < iResidue->nContacts; i++) {
        if(iResidue->contactResidueList[i] == jResidue) {
          //  printf("Contact already knwow between residues indexed %d,%d\n", iResidue->index, jResidue->index);
            return 0;
        }
    }
    //printf("ADDING a new contact between residues indexed %d,%d\n", iResidue->index, jResidue->index);
    iResidue->nContacts++;
    iResidue->contactResidueList = realloc( iResidue->contactResidueList, iResidue->nContacts * sizeof(residue_t*) );
    iResidue->contactResidueList[iResidue->nContacts - 1] = jResidue;
    return 1;
}


double distance(atom_t *iAtom, atom_t *jAtom) {
    return sqrt( (iAtom->x - jAtom->x) * (iAtom->x - jAtom->x) + (iAtom->y - jAtom->y) * (iAtom->y - jAtom->y) + (iAtom->z - jAtom->z) * (iAtom->z - jAtom->z) );
}

void getBoundariesCartesian_DUAL(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, atom_t *minCoor, atom_t *maxCoor) {
    atom_t iMinCoor;
    atom_t iMaxCoor;

    atom_t jMinCoor;
    atom_t jMaxCoor;

    getBoundariesCartesian(iAtomList, iAtom, &iMinCoor, &iMaxCoor);
    getBoundariesCartesian(jAtomList, jAtom, &jMinCoor, &jMaxCoor);

    minCoor->x = iMinCoor.x < jMinCoor.x ? iMinCoor.x : jMinCoor.x;
    minCoor->y = iMinCoor.y < jMinCoor.y ? iMinCoor.y : jMinCoor.y;
    minCoor->z = iMinCoor.z < jMinCoor.z ? iMinCoor.z : jMinCoor.z;

    maxCoor->x = iMaxCoor.x > jMaxCoor.x ? iMaxCoor.x : jMaxCoor.x;
    maxCoor->y = iMaxCoor.y > jMaxCoor.y ? iMaxCoor.y : jMaxCoor.y;
    maxCoor->z = iMaxCoor.z > jMaxCoor.z ? iMaxCoor.z : jMaxCoor.z;


}

void getBoundariesCartesian(atom_t * atomList, int nAtom, atom_t *minCoor, atom_t *maxCoor) {

    minCoor->x = 999999.9;
    minCoor->y = 999999.9;
    minCoor->z = 999999.9;
    maxCoor->x = -999999.9;
    maxCoor->y = -999999.9;
    maxCoor->z = -999999.9;

    for (int i = 0; i < nAtom ; i++) {
        minCoor->x = atomList[i].x < minCoor->x ? atomList[i].x : minCoor->x;
        minCoor->y = atomList[i].y < minCoor->y ? atomList[i].y : minCoor->y;
        minCoor->z = atomList[i].z < minCoor->z ? atomList[i].z : minCoor->z;
        maxCoor->x = atomList[i].x > maxCoor->x ? atomList[i].x : maxCoor->x;
        maxCoor->y = atomList[i].y > maxCoor->y ? atomList[i].y : maxCoor->y;
        maxCoor->z = atomList[i].z > maxCoor->z ? atomList[i].z : maxCoor->z;
    }

}

void cartesianToMesh(atom_t *atom, int *i, int *j, int *k, float step, atom_t minCoor) {
    *i = (int) floor( (atom->x - minCoor.x) / step);
    *j = (int) floor( (atom->y - minCoor.y) / step);
    *k = (int) floor( (atom->z - minCoor.z) / step);
}

void dumpMeshContent(meshContainer_t *meshContainer) {
    mesh_t *mesh = meshContainer->mesh;
    int nCells = meshContainer->nFilled;
    cell_t **filledCells = meshContainer->filledCells;
    cell_t *currCell;

    printf("Mesh dimensions : %d %d %d\n", mesh->iMax, mesh->jMax, mesh->kMax);
    for(int i = 0 ; i < nCells; i++) {
        printf("########ANCHOR CELL\n");
        currCell = filledCells[i];
        dumpCellContent(currCell);
        printf(">>>>>NEIGHBOUR CELLS\n");

    }
}
void dumpCellContent(cell_t *cell) {
    printf("Cells %d %d %d has %d members:\n", cell->i, cell->j, cell->k, cell->memberCount);
    atom_t *atom = cell->members;
    char atomString[81];
    while (atom != NULL) {
        stringifyAtom(atom, atomString);
        printf("\t%s\n", atomString);
        atom = atom->nextCellAtom;
    }
}

meshContainer_t *createMeshContainer(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double step) {
    atom_t minCoor;
    atom_t maxCoor;
    bool dualMode = jAtomList != NULL ? true : false;
    if (dualMode)
        getBoundariesCartesian_DUAL(iAtomList, iAtom, jAtomList, jAtom, &minCoor, &maxCoor);
    else
        getBoundariesCartesian(iAtomList, iAtom, &minCoor, &maxCoor);

#ifdef DEBUG
    printf("Minimal Coordinates %g %g %g\n", minCoor.x, minCoor.y, minCoor.z);
    printf("Maximal Coordinates %g %g %g\n", maxCoor.x, maxCoor.y, maxCoor.z);

#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("Maximal Coordinates %g %g %g\n", maxCoor.x, maxCoor.y, maxCoor.z);
    PySys_WriteStdout("Minimal Coordinates %g %g %g\n", minCoor.x, minCoor.y, minCoor.z);
#endif
#endif

    int iDim = (maxCoor.x - minCoor.x);
    iDim = (iDim + step - 1) / step + 1;
    int jDim = (maxCoor.y - minCoor.y);
    jDim = (jDim + step - 1) / step + 1;
    int kDim = (maxCoor.z - minCoor.z);
    kDim = (kDim + step - 1) / step + 1;

    mesh_t *i_mesh = createMesh(iDim, jDim, kDim);
    cell_t ***grid = i_mesh->grid;

#ifdef DEBUG
    printf("Projecting ... \n");
#endif
    // We store the non-empty cells
    cell_t **filledCells = malloc(i_mesh->n * sizeof(cell_t*));
    int nFilled = 0;

    for (int c = 0 ; c < iAtom ; c++) {
        int i, j, k;
        cartesianToMesh(&iAtomList[c], &i, &j, &k, step, minCoor);
        if (grid[i][j][k].memberCount == 0) {
            /*This cell is non-empty
            register its adress*/
            filledCells[nFilled] = &grid[i][j][k];
            // intialize cell data structure
            grid[i][j][k].members = &iAtomList[c];
            grid[i][j][k].iMembers = &iAtomList[c];
            grid[i][j][k].head = grid[i][j][k].members;
            nFilled++;
        } else {
            grid[i][j][k].head->nextCellAtom = &iAtomList[c];
            grid[i][j][k].head = grid[i][j][k].head->nextCellAtom;
        }
        grid[i][j][k].head->nextCellAtom = NULL;
        grid[i][j][k].memberCount++;
        grid[i][j][k].iMemberCount++;
        iAtomList[c].inCell = &(grid[i][j][k]);
    }

    if (dualMode) {
        for (int c = 0 ; c < jAtom ; c++) {
            int i, j, k;
            cartesianToMesh(&jAtomList[c], &i, &j, &k, step, minCoor);
            if (grid[i][j][k].memberCount == 0) {
                /*This cell is non-empty
                register its adress*/
                filledCells[nFilled] = &grid[i][j][k];
                // intialize cell data structure
                grid[i][j][k].members = &jAtomList[c];
                grid[i][j][k].jMembers = &jAtomList[c];
                grid[i][j][k].head = grid[i][j][k].members;
                nFilled++;
            }
            else if (grid[i][j][k].jMemberCount == 0) { // Already created but not j atom
                grid[i][j][k].jMembers = &jAtomList[c];
                grid[i][j][k].head = grid[i][j][k].jMembers;

            } else {
                grid[i][j][k].head->nextCellAtom = &jAtomList[c];
                grid[i][j][k].head = grid[i][j][k].head->nextCellAtom;
            }
            grid[i][j][k].head->nextCellAtom = NULL;
            grid[i][j][k].memberCount++;
            grid[i][j][k].jMemberCount++;
            jAtomList[c].inCell = &(grid[i][j][k]);
        }
    }

#ifdef DEBUG
    printf("%d atoms projected onto %d cells\n", iAtom + jAtom, nFilled);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("%d atoms projected onto %d cells\n", iAtom + jAtom, nFilled);
#endif
#endif
    meshContainer_t *results = malloc (sizeof(meshContainer_t));
    results->mesh = i_mesh;
    results->filledCells = filledCells;
    results->nFilled = nFilled;
    results->step = step;

    return results;
}

meshContainer_t *createMeshContainer_OLD(atom_t *atomList, int nAtom, double step) {

    atom_t minCoor;
    atom_t maxCoor;
    getBoundariesCartesian(atomList, nAtom, &minCoor, &maxCoor);
#ifdef DEBUG
    printf("Minimal Coordinates %g %g %g\n", minCoor.x, minCoor.y, minCoor.z);
    printf("Maximal Coordinates %g %g %g\n", maxCoor.x, maxCoor.y, maxCoor.z);

#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("Maximal Coordinates %g %g %g\n", maxCoor.x, maxCoor.y, maxCoor.z);
    PySys_WriteStdout("Minimal Coordinates %g %g %g\n", minCoor.x, minCoor.y, minCoor.z);
#endif
#endif
    int iDim = (maxCoor.x - minCoor.x);
    iDim = (iDim + step - 1) / step + 1;
    int jDim = (maxCoor.y - minCoor.y);
    jDim = (jDim + step - 1) / step + 1;
    int kDim = (maxCoor.z - minCoor.z);
    kDim = (kDim + step - 1) / step + 1;

    mesh_t *i_mesh = createMesh(iDim, jDim, kDim);
    cell_t ***grid = i_mesh->grid;

#ifdef DEBUG
    printf("Projecting... \n");
#endif
    // We store the non-empty cells
    cell_t **filledCells = malloc(i_mesh->n * sizeof(cell_t*));
    int nFilled = 0;
// Loop over atoms
    for (int c = 0 ; c < nAtom ; c++) {
        int i, j, k;
        cartesianToMesh(&atomList[c], &i, &j, &k, step, minCoor);
        if (grid[i][j][k].memberCount == 0) {
            /*This cell is non-empty
            register its adress*/
            filledCells[nFilled] = &grid[i][j][k];
            // intialize cell data structure
            grid[i][j][k].members = &atomList[c];
            grid[i][j][k].head = grid[i][j][k].members;
            nFilled++;
        } else {
            grid[i][j][k].head->nextCellAtom = &atomList[c];
            grid[i][j][k].head = grid[i][j][k].head->nextCellAtom;
        }
        grid[i][j][k].head->nextCellAtom = NULL;
        grid[i][j][k].memberCount++;
        atomList[c].inCell = &(grid[i][j][k]);
    }
#ifdef DEBUG
    printf("%d atoms projected onto %d cells\n", nAtom, nFilled);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("%d atoms projected onto %d cells\n", nAtom, nFilled);
#endif
#endif
    meshContainer_t *results = malloc (sizeof(meshContainer_t));
    results->mesh = i_mesh;
    results->filledCells = filledCells;
    results->nFilled = nFilled;
    results->step = step;

    return results;
}

cell_t ** vectorizeMesh(mesh_t *mesh) {
    cell_t ** vectorCells = malloc(mesh->n * sizeof(cell_t*));
    int x = 0;
    for(int i = 0; i < mesh->iMax ; i++) {
        for(int j = 0; j < mesh->jMax ; j++) {
            for(int k = 0; k < mesh->kMax ; k++) {
                vectorCells[x] = &mesh->grid[i][j][k];
                x++;
            }
        }
    }
    return vectorCells;
}


void meshDummy(int a, int b, int c) {

    mesh_t *dum_mesh = createMesh(a, b, c);
    printMesh(dum_mesh);

    cell_t **dummy_filledCells = vectorizeMesh(dum_mesh);

    meshContainer_t *dum_results = malloc (sizeof(meshContainer_t*));
    dum_results->mesh = dum_mesh;
    dum_results->filledCells = dummy_filledCells;
    dum_results->nFilled = dum_mesh->n;

    enumerate(dum_results, -1, NULL, false);

    destroyMeshContainer(dum_results);
    return;
}



// Mov to end of the chain and free residues in a backward motion
residue_t *destroyResidueList(residue_t *residueList) {
#ifdef DEBUG
    printf("Destroying residue list\n");
#endif
    residue_t *curr_residue = residueList;
    while (curr_residue->nextResidueList != NULL) {
        curr_residue = curr_residue->nextResidueList;
    }
    while (curr_residue->prevResidueList != NULL) {
        curr_residue = curr_residue->prevResidueList;
        curr_residue->nextResidueList = destroyResidue(curr_residue->nextResidueList);
    }
    curr_residue = destroyResidue(curr_residue);
#ifdef DEBUG
    printf("Ok\n");
#endif
    return curr_residue;
}

residue_t *destroyResidue(residue_t *residue) {
#ifdef DEBUG
    char residueString[81];
    stringifyResidue(residue, residueString);
    printf ("Destroying residue %s\n", residueString);
#endif
    free(residue->contactResidueList);
    free(residue->resName);
    free(residue->resID);
    free(residue);
    return NULL;
}

atom_t *destroyAtomList(atom_t *atomList, int nAtom) {
    //atom_t *next_atom = NULL;
#ifdef DEBUG
    printf ("Destroying atomList\n");
#endif

    for (int n = 0 ; n < nAtom ; n++)
        destroyAtom(&atomList[n]);
    free(atomList); // Destroy the array

#ifdef DEBUG
    printf("All %d atoms list destroyed\n", nAtom);
#endif
    return NULL;
}

atom_t *destroyAtom(atom_t *atom){
#ifdef DEBUG
    char atomString[81];
    stringifyAtom(atom, atomString);
    printf ("Destroying atom %s\n", atomString);
#endif

    free(atom->resID);
    free(atom->resName);
    free(atom->name);
    //free(atom); No need to free atom structure itself, malloc was operate on an array
#ifdef DEBUG
    printf("Ok\n");
#endif
return NULL;
}

meshContainer_t *destroyMeshContainer(meshContainer_t *container) {
#ifdef DEBUG
    printf ("Destroying Mesh Container\n");
#endif
    free(container->filledCells);
    container->mesh = destroyMesh(container->mesh);
    free(container);
#ifdef DEBUG
    printf("All OK\n");
#endif
    return NULL;
}

mesh_t *destroyMesh(mesh_t *i_mesh) {
#ifdef DEBUG
    printf("Destroying Mesh\n");
#endif
    for (int i = 0 ; i < i_mesh->iMax ; i++) {
        for (int j = 0 ; j < i_mesh->jMax ; j++)
            free(i_mesh->grid[i][j]);
        free(i_mesh->grid[i]);
    }
    free(i_mesh->grid);
    free(i_mesh);
#ifdef DEBUG
    printf("Done\n");
#endif
    return NULL;
}

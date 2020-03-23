#include "molecular_object.h"

void printResidueList(residue_t *residueList) {
    residue_t *curr_residue = residueList;
    printf("PRINTING RESIDUE LIST\n");

    printResidue(curr_residue);
    while(curr_residue->nextResidueList != NULL) {
        curr_residue = curr_residue->nextResidueList;
        printResidue(curr_residue);
    }
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

void printAtomList(atom_t *atomList) {
    printf("PAL\n");
    atom_t *curr_atom = &atomList[0];
    printf("PAL1\n");
    char atomString[81];
    while(curr_atom->nextAtomList != NULL) {
        printf("PAL2\n");
        stringifyAtom(curr_atom, atomString);
        printf("%s\n", atomString);
        curr_atom = curr_atom->nextAtomList;
    }
   
    stringifyAtom(curr_atom, atomString);
    printf("%s\n", atomString);
}

void stringifyAtom(atom_t *atom, char *atomString) {
    sprintf(atomString, "%s %s %s %c %g %g %g", atom->name, atom->resName, atom->resID, atom->chainID, atom->x, atom->y, atom->z);
}
void stringifyResidue(residue_t *residue, char *residueString) {
    sprintf(residueString, "%s %s %c [i=%d, n=%d]", residue->resName, residue->resID, residue->chainID, residue->index, residue->nAtoms);
}

void jsonifyResidue(residue_t *residue, char *jsonString) {
    sprintf(jsonString, "{\"resID\" : \"%s\", \"chainID\":\"%c\"}", residue->resID, residue->chainID);
}

// fixed length JSON array stringrepresentation
void jsonArrayifyAtom(atom_t *atom, char *atomString, bool bCoordinates) {
    if (bCoordinates)
        sprintf(atomString, "[ \"%s\", \"%s\", \"%s\", \"%c\", %g, %g, %g ]", atom->name, atom->resName, atom->resID, atom->chainID, atom->x, atom->y, atom->z);
    else 
        sprintf(atomString, "[ \"%s\", \"%s\", \"%s\", \"%c\" ]", atom->name, atom->resName, atom->resID, atom->chainID);
}
void jsonifyAtomPair(atomPair_t *atomPair, char *jsonString) {
    bool bCoordinates = false;
    char atomA[1024];
    char atomB[1024];
    jsonArrayifyAtom(atomPair->a, atomA, bCoordinates);
    jsonArrayifyAtom(atomPair->b, atomB, bCoordinates);
    sprintf(jsonString, "[ %s, %s, %5.3g }", atomA, atomB, atomPair->dist);
}


residue_t *createResidue(atom_t *atom, int n) {
    printf("NINJA\n");
    residue_t *residue = malloc(sizeof(residue_t));
    residue->nAtoms = n;
    residue->nContacts = 0;
    if (atom->name == NULL)
        printf("OUPS\n");
    printf("NINJA %s\n", atom->name);
    residue->resName = malloc( (strlen(atom->resName) + 1) * sizeof(char) );
    printf("NINJA\n");
    strcpy(residue->resName, atom->resName);
    printf("NINJA\n");
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
    //assert(atomList == NULL);
    printf("ZOUM\n");
    residue_t *residue_root = NULL;
    residue_t *residue_head = NULL;
    atom_t *curr_atom = atomList;//&atomList[0];
    atom_t *prev_atom = NULL;
    int nResidue = 0;
    residue_root = createResidue(curr_atom, nResidue);
    residue_head = residue_root;
printf("ZOUM\n");
    curr_atom->belongsTo = residue_head;
    prev_atom = curr_atom;
    curr_atom = curr_atom->nextAtomList;

    while(curr_atom != NULL) {
        printf("ZOUING\n");
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
    printf("ZIM\n");
    return residue_root;
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

// Move to end of the chain and free residues in a backward motion
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

int CreateAtomListFromPdbContainer(pdbCoordinateContainer_t *pdbCoordinateContainer, atom_t *atomList) {
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    //atom_t *atomList = NULL;
    int nAtom = pdbContainerToArrays(pdbCoordinateContainer, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atomList = readFromArrays(nAtom, x, y, z, chainID, resSeq, resName, atomName);
    
    freeAtomListCreatorBuffers(x, y, z, chainID, resSeq, resName, atomName, nAtom);
    return nAtom;
}

void freeAtomListCreatorBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n) {
    free(x);
    free(y);
    free(z);
    free(chainID);
    for (int i = 0; i < n ; i++) {
        free(resID[i]);
        free(resName[i]);
        free(name[i]);
    }
    free(resID);
    free(resName);
    free(name);
}
// MEMORY ALLOCATION OF atom LIST
atom_t *readFromArrays(int nAtoms, double *x, double *y, double *z, char *chainID, char **resID, char **resName, char **name) {
    atom_t *atomList = malloc(nAtoms * sizeof(atom_t));
    for (int n = 0 ; n < nAtoms ; n++) {
        printf("ZZ%d\n", n);
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

atom_t *legacy_readCoordinates(char *fname, int *_nAtom) {
    atom_t *head = NULL;
    atom_t *old = NULL;
    head = malloc(sizeof(atom_t));
    atom_t *root = head;
    FILE * fp;
    size_t len = 0;
    ssize_t read;
    char * line = NULL;

    fp = fopen(fname, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    char *p;
    char *end;

    int nAtom = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        p = line;
        double buf[3];
        int i = 0;
        for (double f = strtod(p, &end); p != end; f = strtod(p, &end)) {
            buf[i] = f;
            p = end;
            //printf("%f\n", f);
            i++;
        }

        head->x = buf[0];
        head->y = buf[1];
        head->z = buf[2];
        head->nextAtomList = malloc(sizeof(atom_t));
        head = head->nextAtomList;
        head->nextAtomList = NULL;
        nAtom++;
    }
    fclose(fp);
    if (line)
        free(line);
    atom_t *atomArray = malloc (nAtom * sizeof(atom_t));
    head = root;
    int i = 0;
    while (head != NULL) {
        if (i < nAtom) {
            atomArray[i].x = head->x;
            atomArray[i].y = head->y;
            atomArray[i].z = head->z;
        }
        old = head;
        head = head->nextAtomList;
        free(old);
        i++;
    }
    free(head);

    *_nAtom = nAtom;
    return atomArray;
}

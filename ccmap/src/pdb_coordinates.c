#include "pdb_coordinates.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "decoygen.h"

#define MAX_RECORD 200000
/*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/

// Deal w/ atom record with truncated lines by filling last position with whitespaces
pdbCoordinateContainer_t *pdbFileToContainer(char *fileName) {
    pdbCoordinateContainer_t *pdbCoordinateContainer = malloc(sizeof(pdbCoordinateContainer_t));
#ifdef DEBUG
    fprintf(stderr, "Reading pdb file %s\n", fileName);
#endif

    FILE *fp;
    fp = fopen(fileName, "r");



    char lineBuffer[100];

    int atomCount = 0;
    char atomFlag[] = "ATOM  ";

    static char atomRecordBuffer[MAX_RECORD][81];
    char recordType[7];


    int terminalOffset;
    while(fgets(lineBuffer, 100, fp)) {
        memcpy( recordType, &lineBuffer[0], 6 );
        recordType[6] = '\0';
        if (strcmp(recordType, atomFlag) != 0) continue;
        strcpy(atomRecordBuffer[atomCount], lineBuffer);
        atomCount++;
        if(atomCount == 1) {
            for (terminalOffset = 80 ; terminalOffset >= 0 ; terminalOffset--) {
                if(atomRecordBuffer[0][terminalOffset] == '\n') {
                    break;
                }
            }
        }
        for (int i = terminalOffset; i < 81 ; i++) {
            atomRecordBuffer[atomCount - 1][i] = ' ';
        }
        //printf("%s", lineBuffer);
    }
    fclose(fp);
#ifdef DEBUG
    fprintf(stderr, "%d atoms have been read\n", atomCount);
#endif


    pdbCoordinateContainer->atomRecordArray = malloc(atomCount * sizeof(atomRecord_t));
    pdbCoordinateContainer->atomCount = atomCount;
    atomRecord_t *newAtom = NULL;
    for (int i = 0; i < atomCount ; i++) {
        newAtom = &(pdbCoordinateContainer->atomRecordArray[i]);

        createAtomRecord(atomRecordBuffer[i], newAtom);
        //stringifyAtomRecord( newAtom, lineBuffer );


     /*   printf("%s", atomRecordBuffer[i]);
        printf("-->\n%s\n\n", lineBuffer);*/
    }


    return pdbCoordinateContainer;
    /*
*/
}

pdbCoordinateContainer_t *destroyPdbCoordinateContainer(pdbCoordinateContainer_t *pdbCoordinateContainer) {
    free(pdbCoordinateContainer->atomRecordArray);
    free(pdbCoordinateContainer);

    return pdbCoordinateContainer;
}

void createAtomRecord(char *recordString, atomRecord_t *newAtom) {
    char buf[10];

    memcpy( newAtom->recordName, &recordString[0], 6 );
    newAtom->recordName[6] = '\0';

    memcpy( buf, &recordString[6], 5 );
    buf[5] = '\0';
    newAtom->serial = atoi(buf);

    memcpy( newAtom->name, &recordString[12], 4 );
    newAtom->name[4] = '\0';

    newAtom->altLoc = recordString[16];

    memcpy( newAtom->resName, &recordString[17], 3 );
    newAtom->resName[3] = '\0';

    newAtom->chainID = recordString[21];

    memcpy( newAtom->resSeq, &recordString[22], 4 );
    newAtom->resSeq[4] = '\0';

    newAtom->iCode = recordString[26];

    memcpy( buf, &recordString[30], 8 );
    buf[8] = '\0';
    newAtom->x = atof(buf);

    memcpy( buf, &recordString[38], 8 );
    buf[8] = '\0';
    newAtom->y = atof(buf);

    memcpy( buf, &recordString[46], 8 );
    buf[8] = '\0';
    newAtom->z = atof(buf);

    memcpy( buf, &recordString[54], 6 );
    buf[6] = '\0';
    newAtom->occupancy = atof(buf);

    memcpy( buf, &recordString[60], 6 );
    buf[6] = '\0';
    newAtom->tempFactor = atof(buf);

    memcpy( newAtom->element, &recordString[76], 2 );
    newAtom->element[2] = '\0';

    memcpy( newAtom->charge, &recordString[78], 2 );
    newAtom->charge[2] = '\0';

    //return newAtom;
}
/*
    Apply optional euler's angles rotation and or translation to passed pdbCoordinateContainer
    Usually translation is intended to center prior to rotation
*/
void transformPdbCoordinateContainer(pdbCoordinateContainer_t *pdbCoordinateContainer, float *euler, float *translation) {
    atomRecord_t *atom;
    float nX, nY, nZ;

    if (translation != NULL) {
        printf("Applying translation %g, %g, %g\n", translation[0], translation[1], translation[2]);
        for (int i = 0 ; i < pdbCoordinateContainer->atomCount ; i++){
            atom = &pdbCoordinateContainer->atomRecordArray[i];
            float oX = atom->x;
            float oY = atom->y;
            float oZ = atom->z;

            atom->x += translation[0];
            atom->y += translation[1];
            atom->z += translation[2];
            printf("[trans]X: %f -> %f, Y:%f -> %f, Z : %f->%f\n", oX, atom->x, oY, atom->y, oZ, atom->z);


        }
    }

    if (euler != NULL) {
        printf("Applying euler rotations %g, %g, %g\n", euler[0], euler[1], euler[2]);
        for (int i = 0 ; i < pdbCoordinateContainer->atomCount ; i++){
            atom = &pdbCoordinateContainer->atomRecordArray[i];
            rotateAtom(atom->x,  atom->y,  atom->z,
                       &nX    ,  &nY    ,  &nZ    ,
                      euler[0],  euler[1], euler[2]);
            printf("[euler]X: %f -> %f, Y:%f -> %f, Z : %f->%f\n", atom->x, nX, atom->y,nY,atom->z, nZ);

            atom->x = nX;
            atom->y = nY;
            atom->z = nZ;
        }
    }

}


void stringifyAtomRecord(atomRecord_t *atomRecord, char *atomRecordString) {
    //printf("|%s|\n", atomRecord->recordName);
    //atomRecordString = NULL;
    /*char test[81];
    sprintf(test, "%6s",\
            atomRecord->recordName);*/
    sprintf(atomRecordString, "%6s%5d %4s%C%3s %c%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s%-2s",\
            atomRecord->recordName, atomRecord->serial, atomRecord->name, atomRecord->altLoc,\
            atomRecord->resName, atomRecord->chainID, atomRecord->resSeq, atomRecord->iCode,\
            atomRecord->x, atomRecord->y, atomRecord->z, atomRecord->occupancy,\
            atomRecord->tempFactor, atomRecord->element, atomRecord->charge);


}

char *pdbContainerToString(pdbCoordinateContainer_t *pdbCoordinateContainer) {
    char buffer[120];
    stringifyAtomRecord(&pdbCoordinateContainer->atomRecordArray[0], buffer);
    int lineLen = strlen(buffer); // The '\n' while overwrite the place reserved for the '\0', strlen adds +1 for '\0'
    int nLines = pdbCoordinateContainer->atomCount;
    int nChar = lineLen *  nLines + 1; // one for '\0'
    char *pdbString = malloc( nChar * sizeof(char) );
    //memset(pdbString, '\0', sizeof(pdbString));
  //  int ccPos = 0;
    for (int i = 0 ; i < pdbCoordinateContainer->atomCount; i++) {
        stringifyAtomRecord( &pdbCoordinateContainer->atomRecordArray[i], buffer);
        strcpy( &(pdbString[ i * lineLen]), buffer);
        pdbString[ (i + 1) * lineLen  - 1] = '\n';
    }
    pdbString[nChar - 1] = '\0';
    return pdbString;
}

void pdbContainerToFile(pdbCoordinateContainer_t *pdbCoordinateContainer, char *fname, char *mode) {
    FILE *fp;

#ifdef DEBUG
    fprintf(stderr, "Opening %s to write %d atoms\n", fname, pdbCoordinateContainer->atomCount);
#endif

    fp=fopen(fname, mode);
    char buffer[120];
    for (int i = 0 ; i < pdbCoordinateContainer->atomCount; i++) {
        stringifyAtomRecord( &pdbCoordinateContainer->atomRecordArray[i], buffer);
        //printf("->%s", buffer);
        fprintf(fp, "%s\n", buffer);
    }
    fclose(fp);
}

int pdbContainerToArrays(pdbCoordinateContainer_t *pdbCoordinateContainer, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {

    int n = pdbCoordinateContainer->atomCount;

    *x = malloc(n * sizeof(double));
    *y = malloc(n * sizeof(double));
    *z = malloc(n * sizeof(double));
    *chainID = malloc(n * sizeof(char));
    *resID = malloc(n * sizeof(char*));
    *resName = malloc(n * sizeof(char*));
    *name = malloc(n * sizeof(char*));
    for (int i = 0 ; i < n ; i++) {

        (*chainID)[i] = pdbCoordinateContainer->atomRecordArray[i].chainID;
        (*x)[i] = pdbCoordinateContainer->atomRecordArray[i].x;
        (*y)[i] = pdbCoordinateContainer->atomRecordArray[i].y;
        (*z)[i] = pdbCoordinateContainer->atomRecordArray[i].z;
        (*resID)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].resSeq) + 1) * sizeof(char));
        strcpy((*resID)[i], pdbCoordinateContainer->atomRecordArray[i].resSeq);
        (*resName)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].resName) + 1) * sizeof(char));
        strcpy((*resName)[i], pdbCoordinateContainer->atomRecordArray[i].resName);
        (*name)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].name) + 1) * sizeof(char));
        strcpy((*name)[i], pdbCoordinateContainer->atomRecordArray[i].name);
    }
    return n;
}


// We dont add iCode here, see wheter it is added to resSeq ou resName

int legacy_readPdbFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {
    pdbCoordinateContainer_t *pdbCoordinateContainer = pdbFileToContainer(fname);

    int n = pdbCoordinateContainer->atomCount;

    *x = malloc(n * sizeof(double));
    *y = malloc(n * sizeof(double));
    *z = malloc(n * sizeof(double));
    *chainID = malloc(n * sizeof(char));
    *resID = malloc(n * sizeof(char*));
    *resName = malloc(n * sizeof(char*));
    *name = malloc(n * sizeof(char*));
    for (int i = 0 ; i < n ; i++) {

        (*chainID)[i] = pdbCoordinateContainer->atomRecordArray[i].chainID;
        (*x)[i] = pdbCoordinateContainer->atomRecordArray[i].x;
        (*y)[i] = pdbCoordinateContainer->atomRecordArray[i].y;
        (*z)[i] = pdbCoordinateContainer->atomRecordArray[i].z;
        (*resID)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].resSeq) + 1) * sizeof(char));
        strcpy((*resID)[i], pdbCoordinateContainer->atomRecordArray[i].resSeq);
        (*resName)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].resName) + 1) * sizeof(char));
        strcpy((*resName)[i], pdbCoordinateContainer->atomRecordArray[i].resName);
        (*name)[i] = malloc((strlen(pdbCoordinateContainer->atomRecordArray[i].name) + 1) * sizeof(char));
        strcpy((*name)[i], pdbCoordinateContainer->atomRecordArray[i].name);
    }

    destroyPdbCoordinateContainer(pdbCoordinateContainer);

    return n;
}

int legacy_readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {
    double xBuffer[20000];
    double yBuffer[20000];
    double zBuffer[20000];
    char resIDBuffer[20000][6];
    char chainIDBuffer[20000];
    char resNameBuffer[20000][4];
    char nameBuffer[20000][5];

    FILE *fp;
    fp = fopen(fname, "r");
    char string[100];
    int n = 0;
    char * pch;

    while(fgets(string, 100, fp)) {
        pch = strtok (string,",\n");
        int nF = 0;
        while (pch != NULL) {
            if (nF == 0) chainIDBuffer[n] = pch[0];
            if (nF == 1)
                strcpy(resIDBuffer[n], pch);
            if (nF == 2) xBuffer[n] = atof(pch);
            if (nF == 3) yBuffer[n] = atof(pch);
            if (nF == 4) zBuffer[n] = atof(pch);
            if (nF == 5)
                strcpy(resNameBuffer[n], pch);
            if (nF == 6)
                strcpy(nameBuffer[n], pch);

            pch = strtok (NULL, ",\n");
            nF++;
        }
        n++;
    }

    *x = malloc(n * sizeof(double));
    *y = malloc(n * sizeof(double));
    *z = malloc(n * sizeof(double));
    *chainID = malloc(n * sizeof(char));
    *resID = malloc(n * sizeof(char*));
    *resName = malloc(n * sizeof(char*));
    *name = malloc(n * sizeof(char*));
    for (int i = 0 ; i < n ; i++) {
        (*x)[i] = 2.45;
        (*chainID)[i] = chainIDBuffer[i];
        (*x)[i] = xBuffer[i];
        (*y)[i] = yBuffer[i];
        (*z)[i] = zBuffer[i];
        (*resID)[i] = malloc((strlen(resIDBuffer[i]) + 1) * sizeof(char));
        strcpy((*resID)[i], resIDBuffer[i]);
        (*resName)[i] = malloc((strlen(resNameBuffer[i]) + 1) * sizeof(char));
        strcpy((*resName)[i], resNameBuffer[i]);
        (*name)[i] = malloc((strlen(nameBuffer[i]) + 1) * sizeof(char));
        strcpy((*name)[i], nameBuffer[i]);
    }
    fclose(fp);
    return n;
}


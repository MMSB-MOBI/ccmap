#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mesh.h"
#include "encode.h"
#include <unistd.h>
#include "pdb_coordinates.h"
#include <getopt.h>

//#include <stddef.h>
//atom_t *readCoordinates(char *fname, int *_nAtom);

atom_t *readCoordinates(char *fname, int *_nAtom) {
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
        /*printf("Retrieved line of length %zu :\n", read);
        printf("%s", line);*/
        p = line;
        double buf[3];
        int i = 0;
        for (double f = strtod(p, &end); p != end; f = strtod(p, &end)) {
            buf[i] = f;
            /*printf("'%.*s' -> ", (int)(end-p), p);*/
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

void freeBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n) {
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

// We dont add iCode here, see wheter it is added to resSeq ou resName

int readPdbFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {
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



int readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {
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

// PASS Transflag and rotate on the fly
void runSingle( char *fname, float dist, int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***) ) {
    /*  ONE SET OF COORDINATES  */
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    atom_t *atomList = NULL;
    int nAtom = 0;

    printf("Reading coordinates from %s\n", fname);
    nAtom = (*readerFunc)(fname, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atomList = readFromArrays(nAtom, x, y, z, chainID, resSeq, resName, atomName);
    char *ccmap = residueContactMap(atomList, nAtom, dist);

    atomList = destroyAtomList(atomList, nAtom);

    printf("JSON Single ccmap\n%s\n", ccmap);
    freeBuffers(x, y, z, chainID, resSeq, resName, atomName, nAtom);
    free(ccmap);
}

void runDual( char *iFname, char *jFname, float dist,
              int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***)
              ) {
    /*  ONE SET OF COORDINATES  */
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    atom_t *atomList = NULL;
    int nAtom = 0;
    printf("Reading coordinates from %s\n", iFname);
    nAtom = (*readerFunc)(iFname, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atomList = readFromArrays(nAtom, x, y, z, chainID, resSeq, resName, atomName);

    /*  OPTIONAL SECOND SET OF COORDINATES  */
    double *x_other;
    double *y_other;
    double *z_other;
    char *chainID_other;
    char **resSeq_other;
    char **resName_other;
    char **atomName_other;
    atom_t *atomList_other = NULL;
    int nAtom_other = 0;
    unsigned int finalLen=0;
    nAtom_other = (*readerFunc)(jFname, &x_other, &y_other, &z_other, &chainID_other, &resSeq_other, &resName_other, &atomName_other);
    atomList_other = readFromArrays(nAtom_other, x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other);

    int *ccmap = residueContactMap_DUAL(atomList, nAtom, atomList_other, nAtom_other, dist, &finalLen);

// CLEAR
    atomList = destroyAtomList(atomList, nAtom);
    atomList_other = destroyAtomList(atomList_other, nAtom_other);
    freeBuffers(x, y, z, chainID, resSeq, resName, atomName, nAtom);
    freeBuffers(x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other, nAtom_other);
#ifdef DEBUG
    fprintf(stderr, "JSON Dual ccmap\\n");
#endif
    for(int i = 0 ; i < finalLen ; i++)
        printf("%d\n", ccmap[i]);
    free(ccmap);
}


// This function exhausticely list atoms forming contacts
void pdbContainerDualAtomList(float dist, pdbCoordinateContainer_t *pdbCoordinateContainerI,  pdbCoordinateContainer_t *pdbCoordinateContainerJ, char *iTag, char *jTag) {

    fprintf(stderr,"Listing atoms\n");
    /*  ONE SET OF COORDINATES  */
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    atom_t *atomList = NULL;
    int nAtom = 0;

    nAtom = pdbContainerToArrays(pdbCoordinateContainerI, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atomList = readFromArrays(nAtom, x, y, z, chainID, resSeq, resName, atomName);
    int *atomListStatus = malloc(nAtom * sizeof(int));

    /*  OPTIONAL SECOND SET OF COORDINATES  */
    double *x_other;
    double *y_other;
    double *z_other;
    char *chainID_other;
    char **resSeq_other;
    char **resName_other;
    char **atomName_other;
    atom_t *atomList_other = NULL;
    int nAtom_other = 0;
    nAtom_other = pdbContainerToArrays(pdbCoordinateContainerJ, &x_other, &y_other, &z_other, &chainID_other, &resSeq_other, &resName_other, &atomName_other);
    atomList_other = readFromArrays(nAtom_other, x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other);
    int *atomListStatus_other = malloc(nAtom_other * sizeof(int));


    atomListInContact(atomList, nAtom, atomList_other, nAtom_other, dist, atomListStatus, atomListStatus_other);

    char buffer[120];
    printf("%s\n", iTag);
    for (int i = 0; i < nAtom ; i++) {
        if (atomListStatus[i]) {
            stringifyAtomRecord( &pdbCoordinateContainerI->atomRecordArray[i], buffer);
            printf("%s\n", buffer);
        }
    }
    printf("%s\n", jTag);
    for (int j = 0; j < nAtom_other ; j++) {
        if (atomListStatus_other[j]) {
            stringifyAtomRecord( &pdbCoordinateContainerJ->atomRecordArray[j], buffer);
            printf("%s\n", buffer);
        }
    }

    // CLEAR
    atomList = destroyAtomList(atomList, nAtom);
    atomList_other = destroyAtomList(atomList_other, nAtom_other);
    freeBuffers(x, y, z, chainID, resSeq, resName, atomName, nAtom);
    freeBuffers(x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other, nAtom_other);
    free(atomListStatus);
    free(atomListStatus_other);
}


void pdbContainerDualCcmap(float dist, pdbCoordinateContainer_t *pdbCoordinateContainerI,  pdbCoordinateContainer_t *pdbCoordinateContainerJ) {
    /*  ONE SET OF COORDINATES  */
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    atom_t *atomList = NULL;
    int nAtom = 0;

    nAtom = pdbContainerToArrays(pdbCoordinateContainerI, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atomList = readFromArrays(nAtom, x, y, z, chainID, resSeq, resName, atomName);

    /*  OPTIONAL SECOND SET OF COORDINATES  */
    double *x_other;
    double *y_other;
    double *z_other;
    char *chainID_other;
    char **resSeq_other;
    char **resName_other;
    char **atomName_other;
    atom_t *atomList_other = NULL;
    int nAtom_other = 0;
    nAtom_other = pdbContainerToArrays(pdbCoordinateContainerJ, &x_other, &y_other, &z_other, &chainID_other, &resSeq_other, &resName_other, &atomName_other);
    atomList_other = readFromArrays(nAtom_other, x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other);

    unsigned int finalLen=0;
    int *ccmap = residueContactMap_DUAL(atomList, nAtom, atomList_other, nAtom_other, dist, &finalLen);

// CLEAR
    atomList = destroyAtomList(atomList, nAtom);
    atomList_other = destroyAtomList(atomList_other, nAtom_other);
    freeBuffers(x, y, z, chainID, resSeq, resName, atomName, nAtom);
    freeBuffers(x_other, y_other, z_other, chainID_other, resSeq_other, resName_other, atomName_other, nAtom_other);
#ifdef DEBUG
    fprintf(stderr, "JSON Dual ccmap\n");
#endif
    for(int i = 0 ; i < finalLen ; i++)
        printf("%d\n", ccmap[i]);
    free(ccmap);
}

void stringToThreeFloats(char *input, float (*vector)[3]) {
    char *err, *p = input;
    if (input == NULL) return;
    float val;
    int i = 0;
    while (*p) {
        val = strtod(p, &err);
        if (p == err) p++;
        else if ((err == NULL) || (*err == 0)) {
            (*vector)[i] = val;
            i++;
        //    printf("Value: %f\n", val);
            break;
        }
        else {
      //      printf("errValue: %f\n", val);
            p = err + 1;
            (*vector)[i] = val;
            i++;
        }
    }
    //printf("->%g %g %g\n", (*vector)[0], (*vector)[1], (*vector)[2]);

}


void parseTransform(char *eulerString, char *translateString, float (*eulers)[3], float (*translate)[3]){
    if (translateString != NULL){
        stringToThreeFloats(translateString, translate);
        //printf ("Cartesian translation vector %g %g %g\n", (*translate)[0], (*translate)[1], (*translate)[2]);
    } else {
       // translate = NULL;
    }
    if (eulerString != NULL){
        stringToThreeFloats(eulerString, eulers);
        //printf ("Euler's angles values %g %g %g\n", (*eulers)[0], (*eulers)[1], (*eulers)[2]);
    } else {
        //eulers = NULL;
    }

}




/*
    Main program to develop and test C library to manipulate PDB structure in Python 2.7
    --rec ReceptorPdbFile
    --lig LigandPdbFile
    --fmt format (default is PDB)

    --transLigInit Initial x,y,z translation vector to origin
    --eulerLig Euler angle combination to final ligand orientation (intended to be performed w/ ligant centered onto origin)
    --transLig" Final x,y,z translation vector from rotated/centered pose to final ligand pose
    --transRecInit" Initial x,y,z translation vector to origin
    --eulerRecInit", Optional rotation for receptor (intended to be performed w/ ligant centered onto origin)

    --dst Treshold distance to compute contact
*/

int main (int argc, char *argv[]) {

    #ifdef DEBUG
    fprintf(stderr, "*** Debug Mode***\n");
    #endif

    int c;
    char *iFile = NULL;
    char *jFile = NULL;
    char *outFile = NULL;
    char *eulerREC = NULL;
    char *eulerLIG = NULL;
    char *translateREC = NULL;
    char *translateLIG = NULL;
    char *translateLIG2 = NULL;
    int listAtomOnly = 0;
    extern char *optarg;
    extern int optind, optopt, opterr;
    int ligTransflg = 0;
    int recTransflg = 0;
    char *optDist = NULL;
    float eulerAngleREC[3]  = { 0.0, 0.0, 0.0 };
    float translationREC[3] = { 0.0, 0.0, 0.0 };
    float eulerAngleLIG[3]  = { 0.0, 0.0, 0.0 };
    float translationLIG[3] = { 0.0, 0.0, 0.0 };
    float translationLIG2[3] = { 0.0, 0.0, 0.0 };

    pdbCoordinateContainer_t *pdbCoordinateContainerJ = NULL;
    pdbCoordinateContainer_t *pdbCoordinateContainerI = NULL;
//int readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {


    int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***) = NULL;

    readerFunc = &readPdbFile;
    const char    *short_opt = "hla:b:e:t:f:d:w:";
    struct option   long_opt[] =
    {
        {"help",               no_argument, NULL, 'h'},
        {"fmt",          required_argument, NULL, 'f'},
        {"rec",          required_argument, NULL, 'a'},
        {"lig",          required_argument, NULL, 'b'},

        {"transLigInit",          required_argument, NULL, 'x'},
        {"eulerLig",          required_argument, NULL, 'y'},
        {"transLig",          required_argument, NULL, 'z'},

        {"transRecInit",          required_argument, NULL, 'u'},
        {"eulerRecInit",          required_argument, NULL, 'v'},

        {"list",               no_argument, NULL, 'l'},

        {"dist",          required_argument, NULL, 'd'},
        {"dump",          required_argument, NULL, 'w'},
        {NULL,            0,                NULL, 0  }
    };

    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        switch(c) {
            case -1:       /* no more arguments */
            case 0:        /* long options toggles */
            break;

            case 'a':
                iFile = strdup(optarg);
                break;
            case 'b':
                jFile = strdup(optarg);
                break;
            case 'x':
                translateLIG = strdup(optarg);
                break;
            case 'y':
                eulerLIG = strdup(optarg);
                break;
            case 'z':
                translateLIG2 = strdup(optarg);
                break;
            case 'u':
                translateREC = strdup(optarg);
                break;
            case 'v':
                eulerREC = strdup(optarg);
                break;
            case 'l':
                listAtomOnly = 1;
                break;
            case 'd':
                optDist = strdup(optarg);
                break;
            case 'w':
                outFile = strdup(optarg);
                break;
            case 'h':
                printf("Main program to develop and test C library to manipulate PDB structure in Python 2.7\n");
                printf("--rec ReceptorPdbFile\n--lig LigandPdbFile\n--fmt format (default is PDB)\n"),
                printf("--transLigInit Initial x,y,z translation vector to origin\n");
                printf("--eulerLig Euler angle combination to final ligand orientation (intended to be performed w/ ligant centered onto origin)\n");
                printf("--transLig Final x,y,z translation vector from rotated/centered pose to final ligand pose\n");
                printf("--transRecInit Initial x,y,z translation vector to origin\n");
                printf("--eulerRecInit Optional rotation for receptor (intended to be performed w/ ligant centered onto origin\n");
                printf("--dst Treshold distance to compute contact\n");
                printf("--list, just list the atoms forming contacts\n");
                printf("-h, --help                print this help and exit\n");
                printf("\n");
                return(0);
      //      case ':':
            case '?':
                fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
                return(-2);

            default:
                fprintf(stderr, "%s: invalid option -- %c\n", argv[0], c);
                fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
                return(-2);
        }
    }


    if (iFile != NULL)
        pdbCoordinateContainerI = pdbFileToContainer(iFile);

    if (jFile != NULL)
        pdbCoordinateContainerJ = pdbFileToContainer(jFile);

// Centering Receptor, eventual rotation
    if (translateREC != NULL || eulerREC != NULL) {
        printf("Moving receptor to initial position\n");
        parseTransform(eulerREC, translateREC, &eulerAngleREC, &translationREC);
        if (pdbCoordinateContainerI != NULL)
            transformPdbCoordinateContainer(pdbCoordinateContainerI, eulerAngleREC, translationREC);
    }

// Centering Ligand then rotate
     if (translateLIG != NULL || eulerLIG != NULL) {
        printf("Moving ligand to initial position AND rotating to final orientation\n");
        parseTransform(eulerLIG, translateLIG, &eulerAngleLIG, &translationLIG);
        if (pdbCoordinateContainerJ != NULL)
            transformPdbCoordinateContainer(pdbCoordinateContainerJ, eulerAngleLIG, translationLIG);
     }

// Move ligand to final position
    if (translateLIG2 != NULL) {
        printf("Moving ligand to final position\n");
        parseTransform(NULL, translateLIG2, NULL, &translationLIG2);
        if (pdbCoordinateContainerJ != NULL)
            transformPdbCoordinateContainer(pdbCoordinateContainerJ, NULL, translationLIG2);
    }

// No Distance, no ccmap computations
    if (optDist != NULL)
        if(pdbCoordinateContainerI != NULL && pdbCoordinateContainerJ != NULL) {
            if (listAtomOnly)
                pdbContainerDualAtomList(atof(optDist), pdbCoordinateContainerI, pdbCoordinateContainerJ, iFile, jFile);
            else
                pdbContainerDualCcmap(atof(optDist), pdbCoordinateContainerI, pdbCoordinateContainerJ);
    }

// Going out
    if(outFile != NULL) {
        pdbContainerToFile(pdbCoordinateContainerI, outFile, "w");
        if (pdbCoordinateContainerJ != NULL)
            pdbContainerToFile(pdbCoordinateContainerJ, outFile, "a");
    }

#ifdef DEBUG
     fprintf(stderr,"Exiting\n");

     fprintf(stderr, "REC_tr: %g, %g, %g\nLIG_tr: %g, %g, %g\nLIG_euler: %g, %g, %g\nLIG_tr_pose: %g, %g, %g\n", \
            translationREC[0], translationREC[1], translationREC[2], \
            translationLIG[0], translationLIG[1], translationLIG[2], \
            eulerAngleLIG[0],  eulerAngleLIG[1] , eulerAngleLIG[2] , \
            translationLIG2[0], translationLIG2[1], translationLIG2[2]);
#endif

    if (pdbCoordinateContainerI != NULL)
        destroyPdbCoordinateContainer(pdbCoordinateContainerI);
    if (pdbCoordinateContainerJ != NULL)
        destroyPdbCoordinateContainer(pdbCoordinateContainerJ);



    exit(0);
}
/*
    if ( errflg || optDist == NULL || iFile == NULL) {
        fprintf(stderr, "usage: . . . ");
        exit(2);
    }

    double step = atof(optDist);

    if(jFile == NULL)
        runSingle(iFile, step, readerFunc);
    else
        runDual(iFile, jFile, step, readerFunc);
*/

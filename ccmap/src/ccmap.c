#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mesh.h"
#include <unistd.h>
#include "pdb_coordinates.h"
#include <getopt.h>
#include "parameters.h"
#include "atom_mapper.h"

// IMPLEMENT CORRECT CALL TO MESH API

/*void contacMap (pdbCoordinateContainer_t *pdbCoordinateContainerI, pdbCoordinateContainer_t *pdbCoordinateContainerJ,
                    float dist, bool bEncode, bool bAtomic) {
                      
    ccmapView_t *ccmap = atomicContactMap(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, bool bEncoded)
}*/
char *computeCCmap( pdbCoordinateContainer_t *pdbCoordinateContainerI, pdbCoordinateContainer_t *pdbCoordinateContainerJ,
                    float dist, bool bEncode, bool bAtomic, atom_map_t *aMap) {
    ccmapView_t *(*computeMap) (atom_t *, int , atom_t *, int, double, bool, bool) = bAtomic\
                        ? &atomicContactMap
                        : &residueContactMap;
    
    float probeRadius = aMap != NULL ? aMap->probeRadius : 0.0;

    atom_t *iAtomList = NULL;
    atom_t *jAtomList = NULL;
    int iAtom = 0;
    int jAtom = 0;
    iAtomList = CreateAtomListFromPdbContainer(pdbCoordinateContainerI, &iAtom, aMap, probeRadius);
    #ifdef DEBUG
    printAtomList(iAtomList, stderr);
    #endif

    if(pdbCoordinateContainerJ != NULL) 
        jAtomList = CreateAtomListFromPdbContainer(pdbCoordinateContainerJ, &jAtom, aMap, probeRadius);
    printf("Computing %s ccmap with %d pdb records as %s...\n", \
                            bAtomic?"atomic":"residue",\
                            pdbCoordinateContainerJ != NULL?2:1,\
                            bEncode?"integers vector":"JSON");
    ccmapView_t *ccmapView = computeMap(iAtomList, iAtom, jAtomList, jAtom, dist, bEncode, aMap != NULL);
    #ifdef DEBUG
    printf("JSON ccmapView:\n%s\n", ccmapView->asJSON);
    #endif
    if (aMap != NULL) {
        string_t *sasaJson = jsonifySasaResults(ccmapView->sasaResults);
        printf("%s\n", sasaJson->value);
        destroyString(sasaJson);
    }
    char *ccmapAsChar;
    if (bEncode) {
        
        /*
        if(bAtomic){
            char msg[] = "USE encodeLen";
            ccmapAsChar = malloc(1024 * sizeof(char));
            strcpy(ccmapAsChar, msg);
        }
        */
        string_t *encodeVector = createString();
        encodeVector->append(encodeVector, "{ \"type\" : \"encodeVector\", \"data\" : [");
        char buffer[81];
        for (int i = 0; i < ccmapView->encodeLen ; i++) {
            if (i != 0)
                encodeVector->append(encodeVector, ", ");
            sprintf(buffer, "%d",ccmapView->asENCODE[i]);
            encodeVector->append(encodeVector, buffer);
        }
        encodeVector->append(encodeVector, " ]}");
        ccmapAsChar = encodeVector->toChar(encodeVector);
        destroyString(encodeVector);
    } else  {
        ccmapAsChar = malloc( (strlen(ccmapView->asJSON) + 1)  * sizeof(char) );
        strcpy(ccmapAsChar, ccmapView->asJSON);
    }
    destroyCcmapView(ccmapView);
    iAtomList = destroyAtomList(iAtomList, iAtom);
    if (jAtomList != NULL)
        jAtomList = destroyAtomList(jAtomList, jAtom);
    assert(iAtomList == NULL && jAtomList == NULL);
    return ccmapAsChar;
}

void displayPdbRecord(char *fileName){
        pdbCoordinateContainer_t *pdbCoordinateContainer = pdbFileToContainer(fileName);
        char *pdbString = pdbContainerToString(pdbCoordinateContainer);
        printf("%s", pdbString);
        free(pdbString);
        destroyPdbCoordinateContainer(pdbCoordinateContainer);
}

void displayHelp() {
    printf("Main program to develop and test C library to manipulate PDB structure in Python 3.X\n\
    --help   , Print this help\n\
    --rec / a, ReceptorPdbFile\n\
    --lig / b, LigandPdbFile\n\
    --dst / d, Treshold distance to compute contact\n\
    --out / o, Valid File path to write resulting contact map\n\
    --tli / u, (TransLigInit) Initial x,y,z translation vector to origin\n\
    --elg / v,  eulerLig Euler angle combination to final ligand orientation (intended to be performed w/ ligand centered onto origin)\n\
    --tl  / k, transLig Final x,y,z translation vector from rotated/centered pose to final ligand pose\n\
    --tri / m, transRecInit Initial x,y,z translation vector to origin\n\
    --eri / n, eulerRecInit Optional rotation for receptor (intended to be performed w/ ligant centered onto origin)\n\
    --enc / e, ccmap int encoding flag\n\
    --atm / c, atomic ccmap flag\n\
    --dst / d, contact distance threshold\n\
    --dump / w, log dump filepath\n\
    --vrad / r, van der waals atom radi definition file\n\
    --fmt  format (default is PDB)\n"
    );
}
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
    char *optResultFile = NULL;
    extern char *optarg;
    extern int optind, optopt, opterr;
    char *optDist = NULL;
    char *optPrad = NULL;
    float eulerAngleREC[3]  = { 0.0, 0.0, 0.0 };
    float translationREC[3] = { 0.0, 0.0, 0.0 };
    float eulerAngleLIG[3]  = { 0.0, 0.0, 0.0 };
    float translationLIG[3] = { 0.0, 0.0, 0.0 };
    float translationLIG2[3] = { 0.0, 0.0, 0.0 };
    bool bAtomic = false;
    bool bEncode = false;
    bool bASA    = false;
    pdbCoordinateContainer_t *pdbCoordinateContainerJ = NULL;
    pdbCoordinateContainer_t *pdbCoordinateContainerI = NULL;
    atom_map_t               *atomMap             = NULL;                             
    float prad = 1.4;
    char *vradFilePath = NULL;
//int readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name) {

    /*
    int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***) = NULL;
    readerFunc = &readPdbFile;
    */
    const char    *short_opt = "hcesa:b:t:f:d:p:w:o:r:";
    struct option   long_opt[] =
    {
        {"help",               no_argument, NULL, 'h'},
        {"fmt",          required_argument, NULL, 'f'},
        {"rec",          required_argument, NULL, 'a'},
        {"lig",          required_argument, NULL, 'b'},

        {"tli",          required_argument, NULL, 'u'},
        {"elg",          required_argument, NULL, 'v'},
        {"tl",          required_argument, NULL,  'k'},

        {"tri",          required_argument, NULL, 'm'},
        {"eri",          required_argument, NULL, 'n'},
        {"dst",          required_argument, NULL, 'd'},
        {"dump",          required_argument, NULL, 'w'},
        {"enc",               no_argument, NULL, 'e'},
        {"atm",               no_argument, NULL, 'c'},
        {"out",               no_argument, NULL, 'o'},
        {"sasa",               no_argument, NULL, 's'},
        {"prad",          required_argument, NULL, 'p'},
        {"vrad",          required_argument, NULL, 'r'},

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
            case 'd':
                optDist = strdup(optarg);
                break;
            case 'p':
                optPrad = strdup(optarg);
            break;
            case 'w':
                outFile = strdup(optarg);
                break;
            case 'c':
                bAtomic = true;
                break;
            case 'e':
                bEncode = true;
                break;
            case 's':
                bASA = true;
                break;
            case 'o':
                optResultFile = strdup(optarg);
            break;
            case 'r':
                vradFilePath = strdup(optarg);
                break;
            case 'h':
                displayHelp();
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

    if (vradFilePath == NULL && bASA){
        fprintf(stderr, "--vrad and --sasa options are mutually required\n");
        displayHelp();
        exit(0);
    }
    /*
    Parsing mandatory 1st molecule aka "REC"
    Eventually moving it
    */
    if (iFile == NULL) {
        fprintf(stderr, "-a option required\n");
        displayHelp();
        exit(0);
    }
    pdbCoordinateContainerI = pdbFileToContainer(iFile);
    // Centering Receptor, eventual rotation
    if (translateREC != NULL || eulerREC != NULL) {
        printf("Moving receptor to initial position\n");
        parseTransform(eulerREC, translateREC, &eulerAngleREC, &translationREC);
        if (pdbCoordinateContainerI != NULL)
            transformPdbCoordinateContainer(pdbCoordinateContainerI, eulerAngleREC, translationREC);
    }
    /*
    Parsing optional 2nd molecule aka "LIG"
    Eventually moving it
    */
    if (jFile != NULL) {
        pdbCoordinateContainerJ = pdbFileToContainer(jFile);
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
    }

    #ifdef DEBUG 
    fprintf(stderr, "REC_tr: %g, %g, %g\nLIG_tr: %g, %g, %g\nLIG_euler: %g, %g, %g\nLIG_tr_pose: %g, %g, %g\n", \
        translationREC[0], translationREC[1], translationREC[2], \
        translationLIG[0], translationLIG[1], translationLIG[2], \
        eulerAngleLIG[0],  eulerAngleLIG[1] , eulerAngleLIG[2] , \
        translationLIG2[0], translationLIG2[1], translationLIG2[2]);
    #endif

    if (optDist != NULL) {
        double step = atof(optDist);
        char *ccmap = computeCCmap(pdbCoordinateContainerI, pdbCoordinateContainerJ, step, bEncode, bAtomic, NULL);
        FILE *fp = optResultFile != NULL       \
                    ? fopen(optResultFile, "w") \
                    : fopen("ccmap.json", "w");
        fprintf(fp, "%s", ccmap);
        fclose(fp);
        free(ccmap); 
    } else if (bASA) {
        prad = optPrad != NULL ? atof(optPrad) : prad;
        atomMap = readAtomMapperFromFile(vradFilePath, prad);
        //atomMapperPrint(atomMap);
        if (atomMap == NULL)
            exit(1);
        char *ccmap = computeCCmap(pdbCoordinateContainerI, pdbCoordinateContainerJ, 2* (VDW_MAX + prad) , bEncode, bAtomic, atomMap);
        free(ccmap);
        destroyAtomMapper(atomMap);
    } else {
        fprintf(stderr, "No contact distance specified, No cc map computed\n");
    }

    // Dumping transformed coordinates
    if(outFile != NULL) {
        pdbContainerToFile(pdbCoordinateContainerI, outFile, "w");
        if (pdbCoordinateContainerJ != NULL)
            pdbContainerToFile(pdbCoordinateContainerJ, outFile, "a");
    }

 
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




/** DUMPYARD **/
/*

void LEGACY_runDual( char *iFname, char *jFname, float dist,
              int (*readerFunc)(char*, double**, double**, double**, char**, char***, char***, char***)
              ) {
    //  ONE SET OF COORDINATES  
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

    //  OPTIONAL SECOND SET OF COORDINATES 
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
void LEGACY_pdbContainerDualAtomList(float dist, pdbCoordinateContainer_t *pdbCoordinateContainerI,  pdbCoordinateContainer_t *pdbCoordinateContainerJ, char *iTag, char *jTag) {

    fprintf(stderr,"Listing atoms\n");
    //  ONE SET OF COORDINATES  
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

    //  OPTIONAL SECOND SET OF COORDINATES  
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


void LEGACY_pdbContainerDualCcmap(float dist, pdbCoordinateContainer_t *pdbCoordinateContainerI,  pdbCoordinateContainer_t *pdbCoordinateContainerJ) {
    //  ONE SET OF COORDINATES  
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

    //  OPTIONAL SECOND SET OF COORDINATES 
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
*/

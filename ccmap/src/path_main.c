#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mesh.h"
#include <unistd.h>
#include "pdb_coordinates.h"
#include <getopt.h>
#include "pathfinder.h"
#include "mesh_default.h"
#include "miscellaneous.h"

void displayHelp(){
    fprintf(stderr, "path_main -x <ATOM_SELECTOR> -y <ATOM_SELECTOR> -i <PDB_FILE_PATH>\n");
    fprintf(stderr, "options:\n\t -c <LINKER CHAIN_ID>\n\t-s <GRID_STEP>\n");
    fprintf(stderr, "\t-o <OUTPUT_FILE [default:\"structure_path.pdb\"]>\n");
    fprintf(stderr, "\t-t <SEARCH_TYPE [allowed:\"point\" or \"surf\", default:\"point\"]>\n");
}

void main_error(char *msg){
    fprintf(stderr, "%s", msg);
    displayHelp();
    exit(23);
}

// Returns the atom element which matches strings based selector of the form 'LYS:A:29:CA'
atom_t *getAtomFromList(atom_t *list, int i_max, stringList_t *selector) {
    atom_t *curr_atom = NULL;
    bool doStrip = true;
    char chainBuffer[2];
    
    for (int i = 0; i < i_max; i++){
        curr_atom = &list[i];
        
        if (!compareStringToChar(selector->elem[0], curr_atom->resName, doStrip))
            continue;
      
        chainBuffer[0] = curr_atom->chainID;
        chainBuffer[1] = '\0';
      
        if (!compareStringToChar(selector->elem[1], chainBuffer, doStrip))
            continue;
      
        if (!compareStringToChar(selector->elem[2], curr_atom->resID, doStrip))
            continue;
      
        if (!compareStringToChar(selector->elem[3], curr_atom->name, doStrip))
            continue;
    
      return curr_atom;
    }
        
    return NULL;
}

void necklaceThreading(pdbCoordinateContainer_t *pdbContainer, meshContainer_t *meshContainer, path_t *best_walk, char segID) {
     // Necklace threading
        double *pearl_x, *pearl_y, *pearl_z = NULL;
        char *pearl_chainID = NULL;
        char **pearl_resID, **pearl_resName, **pearl_name = NULL;
        createRecordArraysFromPath(best_walk, meshContainer,\
            &pearl_x, &pearl_y, &pearl_z, &pearl_chainID, &pearl_resID, &pearl_resName, &pearl_name,\
            segID);    
        appendArraysToPdbContainer(pdbContainer, best_walk->len, \
            pearl_x, pearl_y, pearl_z, pearl_chainID, pearl_resID, pearl_resName,  pearl_name);
    
#ifdef DEBUG
        char *data = pdbContainerToString(pdbContainer);
        printf("%s\n", data);
        free(data);
#endif
        freeAtomListCreatorBuffers(pearl_x, pearl_y, pearl_z,\
            pearl_chainID, pearl_resID, pearl_resName,  pearl_name, best_walk->len);
        
}

// Debuging utility mapping and backmapping of arbirtray i,j,k coordinates
void inspect(meshContainer_t *meshContainer, int i, int j, int k) {

    printf("Getting center coordinates of cell %d %d %d\n", i, j, k);
    double x, y,z;
    meshToCartesian(meshContainer, i, j, k, &x, &y, &z);
    printf("\t=> %f %f %f\n", x, y, z);
    atom_t *dummy = createBareboneAtom(1, x, y, z, 'A', "   1", "DUM", "  C");
    printf("Projecting again %f %f %f\n", x, y, z);
    cell_t *pjCell = getCellFromAtom(meshContainer, dummy);

    printf("\t=>Lands on %d %d %d\n", pjCell->i, pjCell->j, pjCell->k);
}

void appendVoxelToPdbContainer(pdbCoordinateContainer_t *pdbContainer, meshContainer_t *meshContainer, char vID) {
#ifdef DEBUG
    fprintf(stderr, "appending %d voxels to pdb container\n", meshContainer->nVoxels);
#endif

    cell_t *currCell = NULL;
    int n =  meshContainer->nVoxels;
   
    // Allocate space for array buffers
    double *x_vox      = malloc(n * sizeof(double));
    double *y_vox      = malloc(n * sizeof(double));
    double *z_vox      = malloc(n * sizeof(double));
    char *chainID_vox  = malloc(n * sizeof(char));
    char **resID_vox   = malloc(n * sizeof(char*));
    char **resName_vox = malloc(n * sizeof(char*));
    char **name_vox   = malloc(n * sizeof(char*));
    char baseResNameV[4] = "VOX" ;
    char baseResNameS[4] = "SOX" ;
    char baseName[4]    = "CA ";
    char resSeqBuffer[81];
    int i_v = 0;
// Iterate over mesh elements create a VOX atom per voxel'd cell
    for(int i = 0; i < meshContainer->mesh->iMax ; i++) {
        for(int j = 0; j < meshContainer->mesh->jMax ; j++){ 
            for(int k = 0; k < meshContainer->mesh->kMax ; k++) {
                currCell = &meshContainer->mesh->grid[i][j][k];
                if(!currCell->isInterior && !currCell->isSurface)
                    continue;
                meshToCartesian(meshContainer, currCell->i, currCell->j, currCell->k,\
                                               &x_vox[i_v], &y_vox[i_v], &z_vox[i_v]);
                chainID_vox[i_v] = vID;
                sprintf(resSeqBuffer, "%d", i_v + 1 );
                resID_vox[i_v] = malloc((strlen(resSeqBuffer) + 1) * sizeof(char));
                strcpy(resID_vox[i_v], resSeqBuffer);
                resName_vox[i_v] = malloc(4 * sizeof(char));
                strcpy(resName_vox[i_v], \
                    currCell->isSurface ? baseResNameS : baseResNameV);
                name_vox[i_v] = malloc((strlen(baseName) + 1) * sizeof(char));
                strcpy(name_vox[i_v], baseName);
                i_v++;
    }}}
    
            
    appendArraysToPdbContainer(pdbContainer, n, \
        x_vox, y_vox, z_vox, chainID_vox, resID_vox, resName_vox,  name_vox);
    
#ifdef DEBUG
        fprintf(stderr, "Voxels appended to following pdb record");
        char *data = pdbContainerToString(pdbContainer);
        printf("%s\n", data);
        free(data);
#endif
        freeAtomListCreatorBuffers(x_vox, y_vox, z_vox,\
            chainID_vox, resID_vox, resName_vox, name_vox, n);
}

int main (int argc, char *argv[]) {

    #ifdef DEBUG
    fprintf(stderr, "*** Debug Mode***\n");
    #endif

    int c;
    char necklaceID = 'A';
    char *iFile = NULL;
    char *oFile = NULL;
    extern char *optarg;
    extern int optind, optopt, opterr;
    pdbCoordinateContainer_t *pdbCoordinateContainer = NULL;
    float cellDim = 1.4;
    char *x_selector, *y_selector = NULL; // "RESNAME:RESNUM:ATOM:CHAIN"
    char *optCellDim = NULL;
    stringList_t *x_selectorElem, *y_selectorElem = NULL;
    char bufferLog[81];
    char defSearchType[] = "point";
    char *searchType = NULL;
    char defaultOutFile[] = "structure_path.pdb";
    const char    *short_opt = "hi:x:y:o:s:c:t:d";
    char ERROR_LOG[1024];
    bool dry = false;

    struct option   long_opt[] =
    {
        {"help",              no_argument, NULL, 'h'},
        {"pdb",         required_argument, NULL, 'i'},

        {"from",        required_argument, NULL, 'x'},
        {"to",          required_argument, NULL, 'y'},
        {"out",         required_argument, NULL, 'o'},
        {"sz" ,         required_argument, NULL, 's'},
        {"seg",         required_argument, NULL, 'c'},
        {"type",        required_argument, NULL, 't'},
        {"dry",               no_argument, NULL, 'd'},
        {NULL,            0,               NULL,  0 }
    };

    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        switch(c) {
            case -1:       /* no more arguments */
            case 0:        /* long options toggles */
            break;

            case 'i':               
                iFile = strdup(optarg);
                break;
            case 'o':
                oFile = strdup(optarg);
                break;
            case 's':
                optCellDim = strdup(optarg);
                break;
            case 'c':
                necklaceID = optarg[0];
                break;
            case 'x':
                x_selector = strdup(optarg);
                break;
            case 'y':
                y_selector = strdup(optarg);
                break;
            case 't':
                searchType = strdup(optarg);
                break;
            case 'd':
                dry = true;
                break;
            case 'h':
                displayHelp();
                return(0);
            case '?':
                fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
                return(-2);

            default:
                fprintf(stderr, "%s: invalid option -- %c\n", argv[0], c);
                fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
                return(-2);
        }
    }
    
    if(oFile == NULL)
        oFile = defaultOutFile;
    if(searchType == NULL)
        searchType = defSearchType;
    if(!strcmp(searchType, "point") &&
       !strcmp(searchType, "surf")) {
        sprintf(ERROR_LOG, "Wrong search type \"%s\"\n", searchType);
        main_error(ERROR_LOG);
    }
    cellDim = optCellDim != NULL ? atof(optCellDim) : cellDim;
    if( x_selector == NULL ||
        y_selector == NULL )
        main_error("Please specify atom selectors\n");

    if (iFile == NULL)
        main_error("Please specify a pdb file\n");

    x_selectorElem = splitAndCreateStringList(x_selector, ':');
    y_selectorElem = splitAndCreateStringList(y_selector, ':');
#ifdef DEBUG
    string_t *string;
    for (int i = 0 ; i < x_selectorElem->len ; i++) {
        string = x_selectorElem->elem[i];
        string->dump(string);
    }
    for (int i = 0 ; i < y_selectorElem->len ; i++) {
        string = y_selectorElem->elem[i];
        string->dump(string);
    }
#endif
    pdbCoordinateContainer = pdbFileToContainer(iFile);
    int nAtoms = 0;
    // NULL -> no Fibogrid yet
    atom_t *atomList = CreateAtomListFromPdbContainer(pdbCoordinateContainer, &nAtoms, NULL, 1.4);

    atom_t *xAtom    =  getAtomFromList(atomList, nAtoms, x_selectorElem);
    atom_t *yAtom    =  getAtomFromList(atomList, nAtoms, y_selectorElem);
    
    if (xAtom == NULL) {
        sprintf(bufferLog, "atom selector expression \"%s\" doesnt match any atom record\n", x_selector);
        main_error(bufferLog);
    }
    
    if (yAtom == NULL) {
        sprintf(bufferLog, "atom selector expression \"%s\" doesnt match any atom record\n", y_selector);
        main_error(bufferLog);
    }
    char bufAtom[104];
    stringifyAtom(xAtom, bufAtom);
    printf("Start atom:\t%s\n", bufAtom);
    stringifyAtom(yAtom, bufAtom);
    printf("End atom:\t%s\n", bufAtom);
    
    meshContainer_t *meshContainer = createMeshContainer(atomList, nAtoms, NULL, 0, cellDim);
    
    printf("mesh [%dx%dx%d]created w/ %d filled cells\n",\
         meshContainer->mesh->iMax, meshContainer->mesh->jMax,\
         meshContainer->mesh->kMax, meshContainer->nFilled);
    
    /*
    inspect(meshContainer, 4, 12, 8); // Start cell
    inspect(meshContainer, 5, 12, 8); // Shoud be cell +1
    */
    
   /*if(dry) {
        fprintf(stderr, "DRY RUN:Only building mesh\n");
        exit(1);
    }*/
    //
    path_t *best_walk = searchForPath(meshContainer, searchType,\
        xAtom, yAtom);
    if (strcmp(searchType, "surf") == 0)
        appendVoxelToPdbContainer(pdbCoordinateContainer, meshContainer, 'X');
    
    if (best_walk == NULL) {
        printf("No pathway found connecting specified pair of atoms\n");
    } else {
        printf("###Best pathway\n");
        printf("[%d] %d %d %d\n", best_walk->start->bwfs,\
            best_walk->start->i, best_walk->start->j,\
            best_walk->start->k);
        for(int stp = 0 ; stp < best_walk->len ; stp++)
            printf("[%d] %d %d %d\n", best_walk->cells[stp]->bwfs,\
                best_walk->cells[stp]->i, best_walk->cells[stp]->j,\
                best_walk->cells[stp]->k\
            );
        printf("[%d] %d %d %d\n", best_walk->stop->bwfs,\
            best_walk->stop->i, best_walk->stop->j,\
            best_walk->stop->k);

        necklaceThreading(pdbCoordinateContainer, meshContainer, best_walk, necklaceID);
        best_walk = destroyPath(best_walk);
       
#ifdef DEBUG
        char *pdbAsCharList = pdbContainerToString(pdbCoordinateContainer);
        fprintf(stderr, "Attempting to write following pdb string content to %s\n", oFile);
        fprintf(stderr, "%s", pdbAsCharList);
        free(pdbAsCharList);
#endif
    pdbContainerToFile(pdbCoordinateContainer, oFile, "w");
    }
//Cleanup
        
    destroyStringList(x_selectorElem);
    destroyStringList(y_selectorElem);
    destroyPdbCoordinateContainer(pdbCoordinateContainer);
    destroyAtomList(atomList, nAtoms);
    destroyMeshContainer(meshContainer);
}

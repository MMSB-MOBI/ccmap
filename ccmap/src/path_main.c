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
    fprintf(stderr, "\t---Find the shortest possible path between two atoms ---\n");
    fprintf(stderr, "linky -x <ATOM_SELECTOR> -y <ATOM_SELECTOR> -i <PDB_FILE_PATH>\n");
    fprintf(stderr, "\nNote: ATOM_SELECTOR syntax specifies the following column separated pdb fields surrounded by single quotes:\n");
    fprintf(stderr, "resName:chainID:resNum:name\n\teg: \'LI1:B:301:CA\'\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-o <OUTPUT_FILE [default:\"structure_path.pdb\"]>\n");
    fprintf(stderr, "\t-t <SEARCH_TYPE [allowed:\"point\" or \"surf\", default:\"surf\"]>\n");
    fprintf(stderr, "\t-c <DUMMY_SEGID [default:\'P\']> single char segID of all path atoms in pdb coordinates\n");
    fprintf(stderr, "\t-u <MESH_SIZE [default:0.2]> length of cubic unit cell edges, in A.\n");
    fprintf(stderr, "\t-s <DUM_SPACING [default:3.5]> minimal euclidean distance between path atoms, in A.\n");
    fprintf(stderr, "\t-w <H20_radius [default:1.4]> water probe radius for solvant exclusion surface, in A.\n");
    fprintf(stderr, "\t--force (do not early exit on atom duplicates)\n");
    fprintf(stderr, "\t--pshow (show cells path in output)\n");
    fprintf(stderr, "\t--dry (Only compute the mesh)\n");
    fprintf(stderr, "\t--vshow (show voxels coordinates in output, may cause pdb numbering overflow!!!)\n");
    fprintf(stderr, "\t--zoom (increase Fibonacci grid resolution, trading speed for accuracy)\n");
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

// Reconstruct a thread of CA along polyline w/ decent spacings
void necklaceThreading(pdbCoordinateContainer_t *pdbContainer, \
                       meshContainer_t *meshContainer, path_t *best_walk, char segID,\
                       double euclidStep, int *newCA, double *trailDist) {
                        
        double *pearl_x, *pearl_y, *pearl_z = NULL;
        char *pearl_chainID = NULL;
        char **pearl_resID, **pearl_resName, **pearl_name = NULL;
        int status = createRecordArraysFromPath(best_walk, meshContainer,\
            &pearl_x, &pearl_y, &pearl_z, &pearl_chainID, &pearl_resID, &pearl_resName, &pearl_name,\
            segID, euclidStep, newCA, trailDist);
            if(status == -1) {
                fprintf(stderr, "Recaord Array Fatal Error\n");
                exit(-1); 
            }
        appendArraysToPdbContainer(pdbContainer, *newCA, \
            pearl_x, pearl_y, pearl_z, pearl_chainID, pearl_resID, pearl_resName,  pearl_name);
        
#ifdef DEBUG
        char *data = pdbContainerToString(pdbContainer);
        printf("%s\n", data);
        free(data);
#endif
        freeAtomListCreatorBuffers(pearl_x, pearl_y, pearl_z,\
            pearl_chainID, pearl_resID, pearl_resName,  pearl_name, *newCA);
}

// Debuging utility mapping and backmapping of arbirtray i,j,k coordinates
void inspect(meshContainer_t *meshContainer, int _i, int _j, int _k, atom_t *iAtom) {
    char logAtom[1024];
    int i,j,k;
    cell_t *pjCell = NULL;
    if(iAtom == NULL) {
        printf("\t  --- inspect arbitrary cell start mode---\n");
        i = _i;
        j = _j;
        k = _k;
    } else {
        stringifyAtom(iAtom, logAtom);
        printf("\t --- actual atom start mode ---\nGetting cell coordinates of %s\n", logAtom);
        pjCell = getCellFromAtom(meshContainer, iAtom);
        i = pjCell->i;
        j = pjCell->j;
        k = pjCell->k;
    }
    
    
    printf("\tGetting center coordinates of cell %d %d %d\n", i, j, k);
    double x, y,z;
    atom_t *dummy = NULL;
    meshToCartesian(meshContainer, i, j, k, &x, &y, &z);
    printf("\t=> %f %f %f\n", x, y, z);
    dummy = createBareboneAtom(1, x, y, z, 'A', "   1", "DUM", "  C");
    stringifyAtom(dummy, logAtom);
    printf("\tCreating dummy atom object from this meshToCartesian projection:\n\t%s\n", logAtom);
    printf("\tProjecting back in mesh space\n");
    pjCell = getCellFromAtom(meshContainer, dummy);
    dummy = destroyAtom(dummy);
    printf("\t=>Lands on %d %d %d\n\n", pjCell->i, pjCell->j, pjCell->k);
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
    char necklaceID = 'P';
    char *iFile = NULL;
    char *oFile = NULL;
    char *vMapFile = NULL;
    extern char *optarg;
    extern int optind, optopt, opterr;
    pdbCoordinateContainer_t *pdbCoordinateContainer = NULL;
    float cellDim = 0.2;
    float spacing = 3.5;
    char *x_selector, *y_selector = NULL; // "RESNAME:RESNUM:ATOM:CHAIN"
    char *optCellDim = NULL;
    char *optBeadSpc = NULL;
    stringList_t *x_selectorElem, *y_selectorElem = NULL;
    char bufferLog[81];
    char defSearchType[] = "surf";
    char *searchType = NULL;
    char defaultOutFile[] = "structure_path.pdb";
    const char    *short_opt = "hi:x:y:o:s:u:c:t:w:m:vpfz";
    char ERROR_LOG[1024];
    bool dry = false;
    bool vShow = false;
    bool pShow = false;
    bool bSasaHiRes = false;
    bool noAtomCheck = false;
    char *optH20      = NULL;
    double radiusH20 = 1.4;
    atom_map_t * aMap = NULL;
    struct option   long_opt[] =
    {
        {"help",              no_argument, NULL, 'h'},
        {"pdb",         required_argument, NULL, 'i'},
        {"map",         required_argument, NULL, 'm'},
        {"from",        required_argument, NULL, 'x'},
        {"to",          required_argument, NULL, 'y'},
        {"out",         required_argument, NULL, 'o'},
        {"unit" ,         required_argument, NULL, 'u'},
        {"spacing" ,         required_argument, NULL, 's'},
        {"seg",         required_argument, NULL, 'c'},
        {"type",        required_argument, NULL, 't'},
        {"vshow",             no_argument, NULL, 'v'},
        {"pshow",             no_argument, NULL, 'p'},
        {"dry",               no_argument, NULL, 'd'},
        {"force",             no_argument, NULL, 'f'},
        {"zoom",             no_argument, NULL, 'z'},
        //radiusH20
        {"water",             required_argument, NULL, 'w'},
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
            case 'm':
                vMapFile = strdup(optarg);
                break;
            case 'u':
                optCellDim = strdup(optarg);
                break;
            case 's':
                optBeadSpc = strdup(optarg);
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
            case 'w':
                optH20     = strdup(optarg);
                break;
            case 'd':
                dry = true;
                break;            
            case 'v':
                vShow = true;
                break;
            case 'p':
                pShow = true;
                break;  
            case 'f':
                noAtomCheck = true;
                break;
            case 'z':
                bSasaHiRes = true;
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
    cellDim = optCellDim != NULL ? atof(optCellDim) :  cellDim;
    spacing = optBeadSpc != NULL ? atof(optBeadSpc) :  spacing;
    radiusH20 = optH20   != NULL ? atof(optH20)     : radiusH20;
            
        
    if( x_selector == NULL ||
        y_selector == NULL )
        main_error("Please specify atom selectors\n");

    if (iFile == NULL)
        main_error("Please specify a pdb file\n");
    if (vMapFile == NULL)
        main_error("Please specify a mapper file\n");        
    // radius check its usage and validity w/ later usage ^^
    aMap = readAtomMapperFromFile(vMapFile);
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
    atom_t *atomList = CreateAtomListFromPdbContainer(pdbCoordinateContainer, &nAtoms, aMap, radiusH20, bSasaHiRes);
    printf("Applying H20 probe radius of %g A. to atomic solvant volume exclusion\n", radiusH20);
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
    printf("User atom selection:\n");
    char bufAtom[104];
    stringifyAtom(xAtom, bufAtom);
    printf("\tStart atom:\t\"%s\"\n", bufAtom);
    stringifyAtom(yAtom, bufAtom);
    printf("\tEnd atom :\t\"%s\"\n", bufAtom);
    
    meshContainer_t *meshContainer = createMeshContainer(atomList, nAtoms, NULL, 0, cellDim);
    
    printf("Mesh [%dx%dx%d] created: %d cells contain atoms\n",\
         meshContainer->mesh->iMax, meshContainer->mesh->jMax,\
         meshContainer->mesh->kMax, meshContainer->nFilled);
    
#ifdef DEBUG
    stringifyAtom(xAtom, bufAtom);
    printf("\nForward/Backward mapping assertion w/ user start atom...\n%s\n", bufAtom);
    inspect(meshContainer, 9, 9, 9, xAtom); // Start cell
    printf("\nForward/Backward mapping assertion w/ cell 9/9/9...\n");
    inspect(meshContainer, 9, 9, 9, NULL);
    printf("\nForward/Backward mapping assertion w/ cell 35/75/33...\n");
    inspect(meshContainer, 33, 75, 33, NULL);
    
    /*
    inspect(meshContainer, 4, 12, 8); // Start cell
    inspect(meshContainer, 5, 12, 8); // Shoud be cell +1
    */
#endif
    if(dry) {
        fprintf(stderr, "DRY RUN:Only building mesh\n");
        destroyStringList(x_selectorElem);
        destroyStringList(y_selectorElem);
        destroyPdbCoordinateContainer(pdbCoordinateContainer);
        destroyAtomList(atomList, nAtoms);
        destroyMeshContainer(meshContainer);
        exit(1);
    }
    //
    path_t *best_walk = searchForPath(meshContainer, searchType,\
        xAtom, yAtom, noAtomCheck);
    if (strcmp(searchType, "surf") == 0 && vShow ) {
        fprintf(stderr, "Trying to append voxel to pdb coordinates");
        appendVoxelToPdbContainer(pdbCoordinateContainer, meshContainer, 'X');
    } 
    if (best_walk == NULL) {
        printf("No pathway found connecting specified pair of atoms\n");
    } else {
        printf("Best pathway -- aprox. polyline lengths %g A\n",\
            xAtom->_radiusASA + yAtom->_radiusASA + cellDim * (best_walk->len +1)
        );
        if (pShow) {
            printf("[%d] %d %d %d\n", best_walk->start->bwfs,\
                best_walk->start->i, best_walk->start->j,\
                best_walk->start->k);
            printf("\t---\n");
            for(int stp = 0 ; stp < best_walk->len ; stp++)
                printf("[%d] %d %d %d\n", best_walk->cells[stp]->bwfs,\
                best_walk->cells[stp]->i, best_walk->cells[stp]->j,\
                best_walk->cells[stp]->k\
                );
            printf("\t---\n");
            printf("[%d] %d %d %d\n", best_walk->stop->bwfs,\
                best_walk->stop->i, best_walk->stop->j,\
                best_walk->stop->k);
        }
      
        int nPathBeads;
        double trailingDist; 
        necklaceThreading(pdbCoordinateContainer, meshContainer, best_walk, necklaceID, spacing, &nPathBeads, &trailingDist);
        printf("Approximate linker length %g A.\n", spacing *(nPathBeads-1) + trailingDist + xAtom->_radiusASA + yAtom->_radiusASA);
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

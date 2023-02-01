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
double necklaceThreading(pdbCoordinateContainer_t *pdbContainer, \
                       meshContainer_t *meshContainer, path_t *best_walk, char segID,\
                       double euclidStep, int *newCA, double *trailDist,\
                       atom_t *atomStart, atom_t *atomStop) {
                        
    double *pearl_x, *pearl_y, *pearl_z;
    char *pearl_chainID = NULL;
    char **pearl_resID, **pearl_resName, **pearl_name = NULL;
    int status = createRecordArraysFromPath(best_walk, meshContainer,\
        &pearl_x, &pearl_y, &pearl_z, &pearl_chainID, &pearl_resID, &pearl_resName, &pearl_name,\
        segID, euclidStep, newCA, trailDist, atomStart, atomStop);
        if(status == -1) {
            fprintf(stderr, "Record Array Fatal Error\n");
            exit(-1); 
        }
    // Compute length on this record array
    double len = 0;
    for (int i = 0 ; i < (*newCA) - 1; i++ )
        len += euclideanDistance3(pearl_x[i], pearl_y[i], pearl_z[i], pearl_x[i + 1], pearl_y[i + 1], pearl_z[i + 1]);
    
    appendArraysToPdbContainer(pdbContainer, *newCA, \
        pearl_x, pearl_y, pearl_z, pearl_chainID, pearl_resID, pearl_resName,  pearl_name);
    
#ifdef DEBUG
        char *data = pdbContainerToString(pdbContainer);
        printf("%s\n", data);
        free(data);
#endif
    freeAtomListCreatorBuffers(pearl_x, pearl_y, pearl_z,\
        pearl_chainID, pearl_resID, pearl_resName,  pearl_name, *newCA);
    return len;
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
    int SasaResLvl = 1;
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
                SasaResLvl = 2;
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
    atom_t *atomList = CreateAtomListFromPdbContainer(pdbCoordinateContainer, &nAtoms, aMap, radiusH20, SasaResLvl);
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
    /*stringifyAtom(xAtom, bufAtom);
    printf("\nForward/Backward mapping assertion w/ user start atom...\n%s\n", bufAtom);
    inspect(meshContainer, 9, 9, 9, xAtom); // Start cell
    printf("\nForward/Backward mapping assertion w/ cell 9/9/9...\n");
    inspect(meshContainer, 9, 9, 9, NULL);
    printf("\nForward/Backward mapping assertion w/ cell 35/75/33...\n");
    inspect(meshContainer, 33, 75, 33, NULL);
    */
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
        destroyAtomMapper(aMap);
        exit(1);
    }
    //
    path_t *best_walk = searchForPath(meshContainer, searchType,\
        xAtom, yAtom, noAtomCheck);
    if (strcmp(searchType, "surf") == 0 && vShow ) {
        fprintf(stderr, "Trying to append voxel to pdb coordinates");
        appendVoxelToPdbContainer(pdbCoordinateContainer, meshContainer, 'X', false);
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
        double necklace_len = necklaceThreading(pdbCoordinateContainer, meshContainer, best_walk, \
        necklaceID, spacing, &nPathBeads, &trailingDist,xAtom, yAtom);
        printf("Approximate linker length %g A.\n", necklace_len);
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
    destroyAtomMapper(aMap);        
}

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

int main (int argc, char *argv[]) {

    #ifdef DEBUG
    fprintf(stderr, "*** Debug Mode***\n");
    #endif

    int c;
    char *iFile = NULL;
    char *oFile = NULL;
    extern char *optarg;
    extern int optind, optopt, opterr;
    pdbCoordinateContainer_t *pdbCoordinateContainer = NULL;
    float cell_dim = 1.4;
    char *x_selector, *y_selector = NULL; // "RESNAME:RESNUM:ATOM:CHAIN"
    
    stringList_t *x_selectorElem, *y_selectorElem = NULL;
    char bufferLog[81];
    const char    *short_opt = "hi:x:y:o:s:";
    struct option   long_opt[] =
    {
        {"help",              no_argument, NULL, 'h'},
        {"pdb",         required_argument, NULL, 'i'},

        {"from",        required_argument, NULL, 'x'},
        {"to",          required_argument, NULL, 'y'},
        {"out",         required_argument, NULL, 'o'},
        {"sz" ,         required_argument, NULL, 's'},
       
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
                oFile = strdup(optarg);
                break;

            case 'x':
                x_selector = strdup(optarg);
                break;
            case 'y':
                y_selector = strdup(optarg);
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
    
    if( x_selector == NULL ||
        y_selector == NULL )
        main_error("Please specify atom selectors\n");

    if (iFile == NULL)
        main_error("Please specify a pdb file\n");

    x_selectorElem = splitAndCreateStringList(x_selector, ':');
    y_selectorElem = splitAndCreateStringList(y_selector, ':');
    string_t *string;
    for (int i = 0 ; i < x_selectorElem->len ; i++) {
        string = x_selectorElem->elem[i];
        string->dump(string);
    }
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

    meshContainer_t *meshContainer = createMeshContainer(atomList, nAtoms, NULL, 0, 1.0);
    printf("mesh [%dx%dx%d]created w/ %d filled cells\n",\
         meshContainer->mesh->iMax, meshContainer->mesh->jMax,\
         meshContainer->mesh->kMax, meshContainer->nFilled);
    path_t *best_walk = searchForPath(meshContainer, xAtom, yAtom);
    //computePath(mesh, );
    //destroyStringList(x_selectorElem[0]);


    /*
    Parsing mandatory pdb molecule
    */
   /*
    if (iFile == NULL) {
        fprintf(stderr, "-i option required\n");
        displayHelp();
        exit(0);
    }
    */
   destroyStringList(x_selectorElem);
   destroyStringList(y_selectorElem);
   destroyPdbCoordinateContainer(pdbCoordinateContainer);
   destroyAtomList(atomList, nAtoms);
   destroyMeshContainer(meshContainer);
}

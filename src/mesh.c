#include "mesh.h"


ccmapView_t *atomicContactMap(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, bool bEncoded) {
    bool bAtomic = true;
    ccmapResults_t *ccmapResults = ccmapCore(iAtomList, iAtom, jAtomList, jAtom, ctc_dist, bAtomic);
    ccmapView_t *ccmapView = createCcmapView();
    
    string_t *jsonString = jsonifyAtomPairList(ccmapResults);
    strcpy(ccmapView->asJSON, jsonString->value);
    destroyString(jsonString);
    
    return ccmapView;
}

// ENCODED single residue set currently DISABLED
ccmapView_t *residueContactMap(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, bool bEncoded) {
    bool bAtomic = false;
     fprintf(stderr, "POUET\n");
    ccmapResults_t *ccmapResults = ccmapCore(iAtomList, iAtom, jAtomList, jAtom, ctc_dist, bAtomic);
     fprintf(stderr, "POUET\n");
    ccmapView_t *ccmapView = createCcmapView();
    fprintf(stderr, "POUET\n");
    if (jAtomList != NULL) { 
    // Link the two residues list
        fuseResidueLists(ccmapResults->iResidueList, ccmapResults->jResidueList);
        ccmapResults->fused = true;
    } else {
        assert(!bEncoded); // MUST CHECK FOR encoding one residueList integrity
    }
    if (bEncoded) {
        unsigned int finalLen;
        int jlen = chainLen(ccmapResults->jResidueList);
        int ilen = chainLen(ccmapResults->iResidueList);
        ccmapView->asENCODE = encodeContactMap(ccmapResults->iResidueList, jlen, ilen, &finalLen);
    } else {
        ccmapView->asJSON = jsonifyContactList(ccmapResults->iResidueList);
    }
    destroyCcmapResults(ccmapResults);
    return ccmapView;
}

ccmapView_t *createCcmapView() {
    ccmapView_t *ccmapView = malloc( sizeof(ccmapView_t) );
    ccmapView->asJSON = NULL;
    ccmapView->asENCODE = NULL;
    return ccmapView;
}
ccmapView_t *destroyCcmapView(ccmapView_t *ccmapView) {
    if (ccmapView->asJSON != NULL)
        free(ccmapView->asJSON);
    else
        free(ccmapView->asENCODE);
    free(ccmapView);
    return ccmapView;
}

// IS IT SAFE TO DESTROY jResidue after fusion?
ccmapResults_t *destroyCcmapResults (ccmapResults_t *results){
    results->iResidueList = destroyResidueList(results->iResidueList);
    if (!results->fused)
        results->jResidueList = destroyResidueList(results->jResidueList);
    destroyCellCrawler(results->cellCrawler);
    free(results);
    return results;
}

ccmapResults_t *ccmapCore(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, bool bAtomic) {
    
    residue_t *iResidueList                     = createResidueList(iAtomList);
    fprintf(stderr, "POUM\n");
    residue_t *jResidueList = jAtomList != NULL ? createResidueList(jAtomList) : NULL;
    fprintf(stderr, "POUM\n");
#ifdef DEBUG
#ifdef AS_PYTHON_EXTENSION
        PySys_WriteStdout("Computing residue contact map w/ %.2g Angstrom step\n", ctc_dist);
#endif
    printf("Computing residue contact map w/ %.2g Angstrom step\n", ctc_dist);
#endif
    double step = ctc_dist;
    meshContainer_t *meshContainer = createMeshContainer(iAtomList, iAtom, jAtomList, jAtom, step);
    /* Inspecting atom projection */
    // 101_B_CE1 and 121_1_OE1 cell coordinates ?
    // printResidueCellProjection(" 101", 'B', results, iResidueList);
    //printResidueCellProjection(" 121", 'A', results, jResidueList);
    bool dual = jResidueList != NULL;
    cellCrawler_t cellCrawler = createCellCrawler(bAtomic, dual, ctc_dist);
    meshCrawler(meshContainer, &cellCrawler);
    ccmapResults_t *results = malloc(sizeof(results));
    results->iResidueList = iResidueList;
    results->jResidueList = jResidueList;
    results->cellCrawler = &cellCrawler;
    
    meshContainer = destroyMeshContainer(meshContainer);

    #ifdef DEBUG
        printContactList(iResidueList);
        printf("%s\n", jsonString);
    #ifdef AS_PYTHON_EXTENSION
        PySys_WriteStderr("%s\n", jsonString);
        PySys_WriteStderr("\n" );
    #endif
    #endif

    return results;
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
/*
function
<cstdio>
sprintf
int sprintf ( char * str, const char * format, ... );
Write formatted data to string
Composes a string with the same text that would be printed if format was used on printf, but instead of being printed, the content is stored as a C string in the buffer pointed by str.
The size of the buffer should be large enough to contain the entire resulting string (see snprintf for a safer version).

A terminating null character is automatically appended after the content.
*/
string_t *jsonifyAtomPairList(ccmapResults_t *ccmapResults) {
    fprintf(stderr, "FINALLY\n");
    string_t *jsonString = createString();
    char buffer[1024];
    atomPair_t *atomPair = ccmapResults->cellCrawler->updater->atomContactList;
    uint64_t nAtomPair = ccmapResults->cellCrawler->updater->totalByAtom;
    
    jsonString->append(jsonString, "{ \"type\":\"atomic\", \"data\" : [");
    if (nAtomPair == 0)
        fprintf(stderr, "No atomic contact to jsonify\n");
    for (int i = 0 ; i  < nAtomPair ; i++) {
            jsonifyAtomPair(&(atomPair[i]), buffer);
            jsonString->append(jsonString, buffer);
    }
    jsonString->append(jsonString, "]}");
    return jsonString;
}   

// MESH ITERATOR
// Go through none empty cells
// Get its following cells neighbours
void meshCrawler(meshContainer_t *meshContainer, cellCrawler_t *cellCrawler) {
    mesh_t *mesh = meshContainer->mesh;
    cell_t **cellList = meshContainer->filledCells;
    int nCells = meshContainer->nFilled;

    cell_t ***grid = mesh->grid;
    cell_t *cur_cell;
    int kStart, jStart;
    bool extractBool = cellCrawler->threshold > 0.0;

#ifdef DEBUG
    int nDist = 0;
    int nContacts = 0;
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
                    cellCrawler->enumerator(cellCrawler, cur_cell, &grid[i][j][k]);
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

/* -----------------   DEBUGING FN  ----------------- */
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

// Debugging function to list the cell coordinates of a specified residues projected atoms
void meshDummy(int a, int b, int c) {

    mesh_t *dum_mesh = createMesh(a, b, c);
    printMesh(dum_mesh);

    cell_t **dummy_filledCells = vectorizeMesh(dum_mesh);

    meshContainer_t *dummyMeshContainer = malloc (sizeof(meshContainer_t*));
    dummyMeshContainer->mesh = dum_mesh;
    dummyMeshContainer->filledCells = dummy_filledCells;
    dummyMeshContainer->nFilled = dum_mesh->n;

    cellCrawler_t dummyCellCrawler = createCellCrawler(true, true, -1); //just filled it, TO CHECK TEST
    meshCrawler(dummyMeshContainer, &dummyCellCrawler);

    dummyMeshContainer = destroyMeshContainer(dummyMeshContainer);
    assert(destroyCellCrawler(&dummyCellCrawler) == NULL);
    return;
}

void printResidueCellProjection(char *resID, char chainID, meshContainer_t *meshContainer, residue_t *residueList) {
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


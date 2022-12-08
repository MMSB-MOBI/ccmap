#include "pathfinder.h"
#include <limits.h>

static int offsets[3] = {-1, 0, 1};
static int bestLen = DEFAULT_BWFS;
/*
    meshContainer_t *meshContainer = createMeshContainer(iAtomList, iAtom, jAtomList, jAtom, step);
    pathfinder(meshContainer, cell_start, cell_end);
    
    cellCrawler_t *cellCrawler = createCellCrawler(bAtomic, dual, ctc_dist, bASA);
    meshCrawler(meshContainer, cellCrawler);

*/

path_t *searchForPath(meshContainer_t *meshContainer,\
    char *type, atom_t *atomStart, atom_t *atomStop, bool force) {

#ifdef DEBUG
    fprintf(stderr, "-- Starting search for path -- \"%s\"\n", type);
#endif

    // Basic point search algorithm    
    cell_t *cell_start = atomStart->inCell; 
    cell_t *cell_stop = atomStop->inCell;
    bool (*cellPredicate)(cell_t*) = &pointExplorerPredicate;
    char atomLog[81];
    if(strcmp(type, "surf") == 0 ){ 
        cellPredicate = &surfaceExplorerPredicate;
        printf("Building surfaces w/ mesh unit of %g A. ...\n", \
            meshContainer->step);
        if (!buildSurfaces(meshContainer, force))
            exit(1);
        printf("\tTotal of %d voxels constructed\n", meshContainer->nVoxels);
        
        // TO DO for start and stop
        // Get set of cells that are surface 
        // Identify patches, ie the collection of cells that are connex
        // Set as many start/stop points as there are patches

        printf("Searching for start/stop cells at start/stop atoms surfaces...\n");
        setCells_t *start_cells = getSurfaceCells(atomStart, meshContainer);
        if(start_cells->size == 0) {
            stringifyAtom(atomStart, atomLog);
            fprintf(stderr, "Fatal: no solvant accessible START cell detected for %s\n", atomLog);
            return NULL;
        }
        setCells_t *stop_cells = getSurfaceCells(atomStop, meshContainer);
        if(stop_cells->size == 0) {
            stringifyAtom(atomStop, atomLog);
            fprintf(stderr, "Fatal: no solvant accessible STOP cell detected for %s\n", atomLog);
            return NULL;
        }
        for(int i_stop = 0 ; i_stop < stop_cells->size ; i_stop++)
            stop_cells->cells[i_stop]->isStop = true;
        printf("\tstart/stop surfaces contain %d/%d voxels, one naive cell picking on each...\n",\
            start_cells->size, stop_cells->size);
        cell_start = start_cells->cells[0];
        cell_stop  = stop_cells->cells[0];
        start_cells = destroySetCells(start_cells);
        stop_cells  = destroySetCells(stop_cells);
    }
#ifdef DEBUG
    fprintf(stderr, "Inital best length %d\n", bestLen);
#endif

    cell_start->bwfs = 0;
    printf("Starting from cell (%d,%d, %d) (b=%d)\n", cell_start->i, cell_start->j, cell_start->k,\
    cell_start->memberCount);
    printf("Trying to reachcell (%d,%d, %d) (b=%d)\n", cell_stop->i, cell_stop->j, cell_stop->k,\
    cell_start->memberCount);
  
    exploreCell(meshContainer, cellPredicate,\
                cell_start, 0, cell_start, cell_stop);
    
    if (bestLen >= DEFAULT_BWFS) 
        return NULL;

//#ifdef DEBUG
    printf("\t---Best walk is made of %d moves---\n", bestLen);
//#endif
    
    path_t *path = backtrack(meshContainer, cell_start, cell_stop, type);

    return path;
}

path_t *backtrack(meshContainer_t *meshContainer, cell_t *startCell, cell_t *stopCell, \
                  char *type) {
    path_t *best_path = malloc(sizeof (path_t));

    //printf("STOP CELL BWFS = %d\n", stopCell->bwfs);
    best_path->len = stopCell->bwfs - 1;

    int nbBackSteps = (int)(best_path->len >0?best_path->len:0);
    best_path->cells = malloc( nbBackSteps * sizeof(cell_t*) );
    //fprintf(stderr, "Allocating %d cells path ptr\n", nbBackSteps );
    best_path->start = startCell;
    best_path->stop  = stopCell;
    
    bool (*cellPredicate)(cell_t*) = &pointExplorerPredicate;
    if(strcmp(type, "surf") == 0 )
        cellPredicate = &surfaceExplorerPredicate;

    cell_t *buffer_cell = stopCell;
    for (int i_step = best_path->len - 1 ; i_step >= 0 ; i_step--) {
        buffer_cell = walkBack(buffer_cell, startCell, meshContainer, cellPredicate);    
        best_path->cells[i_step] = buffer_cell;
    }
    return best_path;
}
path_t *destroyPath(path_t *path){
#ifdef DEBUG
    fprintf(stderr, "Destroying path object\n");
#endif
    free(path->cells);
    free(path);
    return NULL;
}

cell_t *walkBack(cell_t *currCell, cell_t *targetCell, meshContainer_t *meshContainer,\
                 bool (*cellPredicate)(cell_t*) ) {
    
#ifdef DEBUG
    fprintf(stderr, "Walking Back through %d %d %d bwfs= %d\n",\
        currCell->i, currCell->j, currCell->k, currCell->bwfs);
#endif

    offsets_t moves[26];
    short int neighbourCount = sortNeighboursByMeshDistanceFace2Face(currCell, \
                                targetCell, meshContainer->mesh, moves);
/*    short int neighbourCount = sortNeighboursByMeshDistance(currCell, \
                                targetCell, meshContainer->mesh, moves);
*/  
    cell_t *closest_cell =  currCell; // seems ok as a dummy initializer
    cell_t *buff_cell    =  NULL;
    for (int iCell = 0 ; iCell < neighbourCount ; iCell++) {
        int i = moves[iCell].i;
        int j = moves[iCell].j;        
        int k = moves[iCell].k;
        buff_cell = &meshContainer->mesh->grid[i][j][k];
        /*
        fprintf(stderr, "Walking Back Assessing %d %d %d bwfs= %d isInterior/isSurface %s/%s\n",\
        buff_cell->i, buff_cell->j, buff_cell->k, buff_cell->bwfs,\
        buff_cell->isInterior?"true":"false", buff_cell->isSurface?"true":"false");
        */
        if( !cellPredicate(buff_cell) )
            continue;
        if (buff_cell->bwfs <  closest_cell->bwfs)
            closest_cell = buff_cell;
    }
#ifdef DEBUG
    fprintf(stderr, "Walking Back TO %d %d %d bwfs= %d\n",\
        closest_cell->i, closest_cell->j, closest_cell->k, closest_cell->bwfs);
#endif
    return closest_cell;
}

bool areSameCells(cell_t *a, cell_t *b) {
    return (a->i == b->i) && (a->j == b->j) && (a->k == b->k);
}

void exploreCell(meshContainer_t *meshContainer, bool (*cellPredicate)(cell_t*),\
    cell_t *currentCell, int nStepFromStart, cell_t *startCell, cell_t *endCell) {
    offsets_t moves[26];
    //20 56 46
    //61 98 77
    // TT : 20 56 49
    int ii = -1;
    int jj = -1;
    int kk = -1;
    //fprintf(stderr, "Hello [%d %d %d](bwfs=%d)\n", currentCell->i, currentCell->j, currentCell->k, currentCell->bwfs);
    if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk) 
        printf("Hello [%d %d %d](bwfs=%d)\n", currentCell->i, currentCell->j, currentCell->k, currentCell->bwfs);
#ifdef DEBUG
    fprintf(stderr, "%d %d %d %d\n", currentCell->i, currentCell->j, currentCell->k, currentCell->bwfs);
#endif
    //printf("%d %d %d %d\n", currentCell->i, currentCell->j, currentCell->k, currentCell->bwfs);
    // Touchdown

    if (areSameCells(currentCell, endCell)) {
#ifdef DEBUG
        fprintf(stderr, "Found destination at depth %d\n", nStepFromStart);
#endif
        bestLen = nStepFromStart;
        currentCell->bwfs = nStepFromStart;
        return;
    }
    // Back to square 1
    if (nStepFromStart > 0 && areSameCells(currentCell, startCell) ) {
        if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk)
            fprintf(stderr, "Back tp sq 1\n");
        return;
    }
    // Cell is blocked
    if (nStepFromStart > 0 && !cellPredicate(currentCell) ) {
        if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk)
            fprintf(stderr, "Cell is blocked\n");
        return; 
    }
    // Exhausted search path
    if (nStepFromStart > 0 && nStepFromStart >= currentCell->bwfs) {
        if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk)
            fprintf(stderr, "Cell is too far\n");
        return;
    }

    if (nStepFromStart > 0)
        if(nStepFromStart < currentCell->bwfs) {
            if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk)
                fprintf(stderr, "Updating bwfs %d %d %d from %d to %d\n",\
                currentCell->i, currentCell->j, currentCell->k,  currentCell->bwfs, nStepFromStart);
            currentCell->bwfs = nStepFromStart;
        }
    // About to exhaust search path
    if (manh_dist(currentCell,  endCell) + currentCell->bwfs > bestLen) {
        if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk) {
            fprintf(stderr, "!! Leaving cell %d [cd(%d %d %d, %d %d %d)] + %d > %d\n", manh_dist(currentCell,  endCell),\
             currentCell->i, currentCell->j, currentCell->k, \
             endCell->i, endCell->j, endCell->k,\
             currentCell->bwfs, bestLen);
        }
#ifdef DEBUG     
        fprintf(stderr, "%d %d %d to far (%d)\n", \
        currentCell->i, currentCell->j, currentCell->k,\
        manh_dist(currentCell,  endCell) + currentCell->bwfs);
#endif
        return;
    }
        
    // Keep on searching
    short int neighbourCount = sortNeighboursByMeshDistanceFace2Face(currentCell, \
                                    endCell, meshContainer->mesh, moves);
    //short int neighbourCount = sortNeighboursByMeshDistance(currentCell, endCell, meshContainer->mesh, moves);
    
    if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk)
        printf(">>[%d %d %d](bwfs=%d) cell feature %d neighbours\n", \
            currentCell->i, currentCell->j, currentCell->k, currentCell->bwfs, neighbourCount);
    for (int iCell = 0 ; iCell < neighbourCount ; iCell++) {
        //61 98 77
        if(currentCell->i == ii && currentCell->j == jj && currentCell->k == kk) {
            printf("[%d %d %d](bwfs=%d) About to move to %d %d %d\n", \
            currentCell->i, currentCell->j, currentCell->k, currentCell->bwfs,\
            moves[iCell].i, moves[iCell].j, moves[iCell].k);
        }

        int i = moves[iCell].i;
        int j = moves[iCell].j;        
        int k = moves[iCell].k;
        exploreCell(meshContainer, cellPredicate,\
                    &meshContainer->mesh->grid[i][j][k],\
                    nStepFromStart + 1, startCell, endCell);
    }

}

bool pointExplorerPredicate(cell_t *cell){
    return cell->memberCount == 0;
}

// sort neigbours by their mesh distance from target destination
// We represent neighbours by their increment
// (+1, +1 , +1) ... (-1, -1, -1)

// Here we only allow face to face transitions to avoid artifacts
short int sortNeighboursByMeshDistance(cell_t *currentCell, cell_t *goal, mesh_t *mesh, offsets_t moves[]){
#ifdef DEBUG
    fprintf(stderr, "\n\t\tSorting neighbourhood\n");
#endif
    int n = 0;
    int new_i, new_j, new_k = 0;
    int i,j,k = 0;
    for (int _i = 0; _i < 3 ; _i++) {
        for (int _j = 0; _j < 3 ; _j++) {
            for (int _k = 0; _k < 3 ; _k++) {
                i = offsets[_i];
                j = offsets[_j];
                k = offsets[_k];

                if (i == 0 && j == 0 && k== 0)
                    continue;
                new_i = currentCell->i + i;
                new_j = currentCell->j + j;
                new_k = currentCell->k + k;
                
                //fprintf(stderr, "yo %d\n", n);
                if(new_i < 0 || new_j < 0 || new_k < 0)
                    continue;
                if(new_i >= mesh->iMax || new_j >= mesh->jMax || new_k >= mesh->kMax)
                    continue;
                    //fprintf(stderr, "??? %d, %d, %d\n", new_i, new_j, new_k);
                /*
                if (mesh->grid[new_i][new_j][new_k].memberCount > 0)
                    continue;
                */
#ifdef DEBUG
                fprintf(stderr, "Registering %d %d %d at pos %d\n", new_i, new_j, new_k, n);
#endif
                moves[n].i = new_i;
                moves[n].j = new_j;
                moves[n].k = new_k;
                moves[n].abs_dist = sqrt( ( new_i - goal->i ) * ( new_i - goal->i ) \
                                        + ( new_j - goal->j ) * ( new_j - goal->j ) \
                                        + ( new_k - goal->k ) * ( new_k - goal->k ) \
                                        );
                n++;
            }
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Total possibile moves in current cell(%d %d %d) is %d\n",\
                currentCell->i,  currentCell->j,  currentCell->k, n);
    for (int x = 0 ; x < n ; x++)
        fprintf(stderr, "[%d] %d %d %d = %f\n", x, moves[x].i, moves[x].j, moves[x].k, moves[x].abs_dist);
    fprintf(stderr, "QSORT\n");
#endif
    qsort(moves, n, sizeof(offsets_t), cmpOffsetfunc);
#ifdef DEBUG
    for (int x = 0 ; x < n ; x++)
        fprintf(stderr, "[%d] %d %d %d = %f\n", x, moves[x].i, moves[x].j, moves[x].k, moves[x].abs_dist);
    fprintf(stderr, "#Sorting ends\n");
#endif
    return n;
}

short int sortNeighboursByMeshDistanceFace2Face(cell_t *currentCell, cell_t *goal, mesh_t *mesh, offsets_t moves[]){
#ifdef DEBUG
    fprintf(stderr, "\n\t\tSorting neighbourhood face2ace\n");
#endif
    int n = 0;
    int new_i, new_j, new_k = 0;
    int i,j,k = 0;
    for (int _i = 0; _i < 3 ; _i++) {
        for (int _j = 0; _j < 3 ; _j++) {
            for (int _k = 0; _k < 3 ; _k++) {
                i = offsets[_i];
                j = offsets[_j];
                k = offsets[_k];
                // face2face transitions
                if ( ! ( (i == 0 && j == 0) ||\
                         (i == 0 && k == 0) ||\
                         (j == 0 && k == 0) \
                    ) )
                    continue;

                if (i == 0 && j == 0 && k== 0)
                    continue;
                new_i = currentCell->i + i;
                new_j = currentCell->j + j;
                new_k = currentCell->k + k;
                
                //fprintf(stderr, "yo %d\n", n);
                if(new_i < 0 || new_j < 0 || new_k < 0)
                    continue;
                if(new_i >= mesh->iMax || new_j >= mesh->jMax || new_k >= mesh->kMax)
                    continue;
                    //fprintf(stderr, "??? %d, %d, %d\n", new_i, new_j, new_k);
                /*
                if (mesh->grid[new_i][new_j][new_k].memberCount > 0)
                    continue;
                */
#ifdef DEBUG
                fprintf(stderr, "Registering %d %d %d at pos %d\n", new_i, new_j, new_k, n);
#endif
                moves[n].i = new_i;
                moves[n].j = new_j;
                moves[n].k = new_k;
                moves[n].abs_dist = sqrt( ( new_i - goal->i ) * ( new_i - goal->i ) \
                                        + ( new_j - goal->j ) * ( new_j - goal->j ) \
                                        + ( new_k - goal->k ) * ( new_k - goal->k ) \
                                        );
                n++;
            }
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Total possibile moves in current cell(%d %d %d) is %d\n",\
                currentCell->i,  currentCell->j,  currentCell->k, n);
    for (int x = 0 ; x < n ; x++)
        fprintf(stderr, "[%d] %d %d %d = %f\n", x, moves[x].i, moves[x].j, moves[x].k, moves[x].abs_dist);
    fprintf(stderr, "QSORT\n");
#endif
    qsort(moves, n, sizeof(offsets_t), cmpOffsetfunc);
#ifdef DEBUG
    for (int x = 0 ; x < n ; x++)
        fprintf(stderr, "[%d] %d %d %d = %f\n", x, moves[x].i, moves[x].j, moves[x].k, moves[x].abs_dist);
    fprintf(stderr, "#Sorting ends\n");
#endif
    return n;
}


// check order
int cmpOffsetfunc (const void * a, const void * b) {
    offsets_t *oA = (offsets_t *)a;
    offsets_t *oB = (offsets_t *)b;
    
    if (oA->abs_dist < oB->abs_dist)
        return -1;
    if (oA->abs_dist > oB->abs_dist)
        return 1;
    return 0;
}
// Generate a necklace of dummy atoms along the path
// Lame code duplication from pdb_coordinates.c: pdbContainerToArrays
// We should clean mem on realloc error exit
int createRecordArraysFromPath(path_t *self, meshContainer_t *meshContainer, double **x, double **y, double **z,\
                                char **chainID, char ***resID, char ***resName, char ***name,\
                                char uID, double spacing) {
    // Theoritical max number of created particules is the length of the path
    int maxNewAtomCount = self->len > 0 ? self->len  : 0;
    
    *x = malloc(maxNewAtomCount * sizeof(double));
    *y = malloc(maxNewAtomCount * sizeof(double));
    *z = malloc(maxNewAtomCount * sizeof(double));
    *chainID = malloc(maxNewAtomCount * sizeof(char));
    *resID = malloc(maxNewAtomCount * sizeof(char*));
    *resName = malloc(maxNewAtomCount * sizeof(char*));
    *name = malloc(maxNewAtomCount * sizeof(char*));

    double *double_swaper    = NULL;
    char *char_swaper        = NULL;
    char **char_array_swaper = NULL;

    double x_prime, y_prime, z_prime;
    char resSeqBuffer[81];
    char baseResName[] = "DUM" ;
    char baseName[] = "CA ";


#ifdef DEBUG
        fprintf(stderr, "Create Necklace atomic parameter arrays over %d cells w/ %g spacing\n"\
    , maxNewAtomCount, spacing);
#endif

    int newAtomCount = 0;
    
    float last_x, last_y, last_z;
    // We consider 1st path cell as safe for particule creation   
    for (int iElem = 0 ; iElem < maxNewAtomCount ; iElem++) {
        meshToCartesian(meshContainer, self->cells[iElem]->i, self->cells[iElem]->j, self->cells[iElem]->k,\
                                       &x_prime             , &y_prime             , &z_prime);
        if (iElem > 0)
            if (euclideanDistance3(x_prime, y_prime, z_prime, last_x, last_y, last_z) < spacing)
                continue;
#ifdef DEBUG        
        if (iElem > 0)
            printf("Threading based on %g distance\n", euclideanDistance3(x_prime, y_prime, z_prime, last_x, last_y, last_z));
#endif
        (*chainID)[newAtomCount] = uID;
        (*x)[newAtomCount] = x_prime;
        (*y)[newAtomCount] = y_prime;
        (*z)[newAtomCount] = z_prime;
        sprintf(resSeqBuffer, "%d", newAtomCount + 1 );
        (*resID)[newAtomCount] = malloc((strlen(resSeqBuffer) + 1) * sizeof(char));
        strcpy((*resID)[newAtomCount], resSeqBuffer);
        (*resName)[newAtomCount] = malloc((strlen(baseResName) + 1) * sizeof(char));
        strcpy((*resName)[newAtomCount], baseResName);
        (*name)[newAtomCount] = malloc((strlen(baseName) + 1) * sizeof(char));
        strcpy((*name)[newAtomCount], baseName);
#ifdef DEBUG
        fprintf(stderr, "%f %f %f \'%s\' \'%s\' \'%c\' \'%s\'\n",\
            (*x)[newAtomCount], (*y)[newAtomCount], (*z)[newAtomCount],\
            (*resName)[newAtomCount], (*resID)[newAtomCount], \
            (*chainID)[newAtomCount], (*name)[newAtomCount] );
#endif
        last_x = x_prime;
        last_y = y_prime;
        last_z = z_prime;
        newAtomCount++;
    }

    double_swaper = (double *) realloc(*x, sizeof(double) * newAtomCount);
    if(double_swaper == NULL) {
        reallocErrorLog(maxNewAtomCount, newAtomCount, "x coordinates");
        return -1;
    }
    *x = double_swaper;

    double_swaper = (double *) realloc(*y, sizeof(double) * newAtomCount);
    if(double_swaper == NULL) {
        reallocErrorLog(maxNewAtomCount, newAtomCount, "y coordinates");
        return -1;
    }
    *y = double_swaper;

    double_swaper = (double *) realloc(*z, sizeof(double) * newAtomCount);
    if(double_swaper == NULL) {
       reallocErrorLog(maxNewAtomCount, newAtomCount, "z coordinates");
        return -1;
    }
    *z = double_swaper;

    char_swaper = (char *) realloc(*chainID, sizeof(char) * newAtomCount);
    if(char_swaper == NULL) {
        reallocErrorLog(maxNewAtomCount, newAtomCount, "segID");
        return -1;
    }
    *chainID = char_swaper;

    char_array_swaper = (char **) realloc(*resID, sizeof(char*) * newAtomCount);
    if(char_array_swaper == NULL) {
        reallocErrorLog(maxNewAtomCount, newAtomCount, "resID");
        return -1;
    }
    *resID = char_array_swaper;

    char_array_swaper = (char **) realloc(*resName, sizeof(char*) * newAtomCount);
    if(char_array_swaper == NULL) {
        reallocErrorLog(maxNewAtomCount, newAtomCount, "resName");
        return -1;
    }
    *resName = char_array_swaper;

    char_array_swaper = (char **) realloc(*name, sizeof(char*) * newAtomCount);
    if(char_array_swaper == NULL) {
        reallocErrorLog(maxNewAtomCount, newAtomCount, "name");
        return -1;
    }
    *name = char_array_swaper;
    printf("Trailing space equals %gA\n", euclideanDistance3( last_x,  last_y, last_z,\
                                                            x_prime, y_prime, z_prime));
    printf("Threading of %d atoms w/ %gA spacing completed\n", newAtomCount, spacing);
    
    return newAtomCount;
}

void reallocErrorLog(int a, int b, char type[]) {
    fprintf(stderr, "Could not reallocation memory for atomRecord %s list from %d to%d\n",\
            type, a, b);
}

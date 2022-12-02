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
    char *type, atom_t *atomStart, atom_t *atomStop) {
    //static unsigned int bestLen = 999999; 
    cell_t *cell_start = atomStart->inCell; 
    cell_t *cell_stop = atomStop->inCell;
    bool (*cellPredicate)(cell_t*) = &pointExplorerPredicate;

    if(strcmp(type, "surf") ){ 
        cellPredicate = &surfaceExplorerPredicate;
        printf("Building surfaces w/ mesh unit= %g ...\n", \
            meshContainer->step);
        buildSurfaces(meshContainer);
        printf("Total of %d voxels constructed\n", meshContainer->nVoxels);
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

#ifdef DEBUG
    fprintf(stderr, "Best length is %d\n", bestLen);
#endif

    return backtrack(meshContainer, cell_start, cell_stop );
}

path_t *backtrack(meshContainer_t *meshContainer, cell_t *startCell, cell_t *stopCell ) {
    path_t *best_path = malloc(sizeof (path_t));
    best_path->len = stopCell->bwfs - 1;
    best_path->cells = malloc( best_path->len * sizeof(cell_t*) );
    best_path->start = startCell;
    best_path->stop  = stopCell;

    cell_t *buffer_cell = stopCell;
    for (int i_step = best_path->len - 1 ; i_step >= 0 ; i_step--) {
       buffer_cell = walkBack(buffer_cell, startCell, meshContainer);
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

cell_t *walkBack(cell_t *currCell, cell_t *targetCell, meshContainer_t *meshContainer) {

#ifdef DEBUG
    fprintf(stderr, "Walking Back through %d %d %d bwfs= %d\n",\
        currCell->i, currCell->j, currCell->k, currCell->bwfs);
#endif

    offsets_t moves[26];
    short int neighbourCount = sortNeighboursByMeshDistance(currCell, \
                                targetCell, meshContainer->mesh, moves);
    
    cell_t *closest_cell =  currCell; // seems ok as a dummy initializer
    cell_t *buff_cell    =  NULL;
    for (int iCell = 0 ; iCell < neighbourCount ; iCell++) {
        int i = moves[iCell].i;
        int j = moves[iCell].j;        
        int k = moves[iCell].k;
        buff_cell = &meshContainer->mesh->grid[i][j][k];
        if (buff_cell->bwfs <  closest_cell->bwfs)
            closest_cell = buff_cell;
    }
    return closest_cell;
}

bool areSameCells(cell_t *a, cell_t *b) {
    return (a->i == b->i) && (a->j == b->j) && (a->k == b->k);
}

void exploreCell(meshContainer_t *meshContainer, bool (*cellPredicate)(cell_t*),\
    cell_t *currentCell, int nStepFromStart, cell_t *startCell, cell_t *endCell) {
    offsets_t moves[26];
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
    if (nStepFromStart > 0 && areSameCells(currentCell, startCell) )
        return;
    // Cell is blocked
    if (nStepFromStart > 0 && cellPredicate(currentCell) )
        return; 
    // Exhausted search path
    if (nStepFromStart > 0 && nStepFromStart >= currentCell->bwfs)
        return;

    if (nStepFromStart > 0)
        currentCell->bwfs = nStepFromStart;
    
    // About to exhaust search path
    if (c_dist(currentCell,  endCell) + currentCell->bwfs > bestLen) {
#ifdef DEBUG     
        fprintf(stderr, "%d %d %d to far (%f)\n", \
        currentCell->i, currentCell->j, currentCell->k,\
        c_dist(currentCell,  endCell) + currentCell->bwfs);
#endif
        return;
    }
    // Keep on searching
    short int neighbourCount = sortNeighboursByMeshDistance(currentCell, endCell, meshContainer->mesh, moves);
    for (int iCell = 0 ; iCell < neighbourCount ; iCell++) {
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
void createRecordArraysFromPath(path_t *self, meshContainer_t *meshContainer, double **x, double **y, double **z,\
                                char **chainID, char ***resID, char ***resName, char ***name,\
                                char uID) {
    int n = self->len;
    *x = malloc(n * sizeof(double));
    *y = malloc(n * sizeof(double));
    *z = malloc(n * sizeof(double));
    *chainID = malloc(n * sizeof(char));
    *resID = malloc(n * sizeof(char*));
    *resName = malloc(n * sizeof(char*));
    *name = malloc(n * sizeof(char*));
    double x_prime, y_prime, z_prime;
    char resSeqBuffer[81];
    char baseResName[] = "DUM" ;
    char baseName[] = "CA ";

#ifdef DEBUG
    fprintf(stderr, "Create Following Necklace arrays [%d elem]\n", n);
#endif

    for (int iElem = 0 ; iElem < n ; iElem++) {
        meshToCartesian(meshContainer, self->cells[iElem]->i, self->cells[iElem]->j, self->cells[iElem]->k,\
                                       &x_prime             , &y_prime             , &z_prime);
        (*chainID)[iElem] = uID;
        (*x)[iElem] = x_prime;
        (*y)[iElem] = y_prime;
        (*z)[iElem] = z_prime;
        sprintf(resSeqBuffer, "%d", iElem + 1 );
        (*resID)[iElem] = malloc((strlen(resSeqBuffer) + 1) * sizeof(char));
        strcpy((*resID)[iElem], resSeqBuffer);
        (*resName)[iElem] = malloc((strlen(baseResName) + 1) * sizeof(char));
        strcpy((*resName)[iElem], baseResName);
        (*name)[iElem] = malloc((strlen(baseName) + 1) * sizeof(char));
        strcpy((*name)[iElem], baseName);
#ifdef DEBUG
        fprintf(stderr, "%f %f %f \'%s\' \'%s\' \'%c\' \'%s\'\n",\
            (*x)[iElem], (*y)[iElem], (*z)[iElem],\
            (*resName)[iElem], (*resID)[iElem], \
            (*chainID)[iElem], (*name)[iElem] );
#endif
    }
}

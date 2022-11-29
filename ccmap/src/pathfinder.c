#include "pathfinder.h"
#include <limits.h>

static int offsets[3] = {-1, 0, 1};
static int bestLen = 999999;
/*
    meshContainer_t *meshContainer = createMeshContainer(iAtomList, iAtom, jAtomList, jAtom, step);
    pathfinder(meshContainer, cell_start, cell_end);
    
    cellCrawler_t *cellCrawler = createCellCrawler(bAtomic, dual, ctc_dist, bASA);
    meshCrawler(meshContainer, cellCrawler);

*/

path_t *searchForPath(meshContainer_t *meshContainer, atom_t *atomStart, atom_t *atomStop) {
    //static unsigned int bestLen = 999999; 
    cell_t *cell_start = atomStart->inCell; 
    cell_t *cell_stop = atomStop->inCell;


    printf("Inital best length %d\n", bestLen);

    cell_start->bwfs = 0;
    printf("Starting from cell (%d,%d, %d) (b=%d)\n", cell_start->i, cell_start->j, cell_start->k,\
    cell_start->memberCount);
    printf("Trying to reachcell (%d,%d, %d) (b=%d)\n", cell_stop->i, cell_stop->j, cell_stop->k,\
    cell_start->memberCount);
    /*
    offsets_t moves[26];
    short int neighbourCount = sortNeighboursByMeshDistance(cell_start, cell_stop, meshContainer->mesh, moves);

    for (int iCell = 0 ; iCell < neighbourCount ; iCell++) {
        int i = moves[iCell].i;
        int j = moves[iCell].j;        
        int k = moves[iCell].k;
        //printf("Browsable at  %d %d %d at c_dist %f\n", i, j, k, moves[iCell].abs_dist);
        exploreCell(meshContainer, &meshContainer->mesh->grid[i][j][k], 0, cell_start, cell_stop);
    }*/
    exploreCell(meshContainer, cell_start, 0, cell_start, cell_stop);
    printf("Best length is %d\n", bestLen);

    return backtrack(meshContainer, cell_start, cell_stop );
}

path_t *backtrack(meshContainer_t *meshContainer, cell_t *startCell, cell_t *stopCell ) {
    path_t *best_path = malloc(sizeof (path_t));
    return best_path;
}

bool areSameCells(cell_t *a, cell_t *b) {
    return (a->i == b->i) && (a->j == b->j) && (a->k == b->k);
}

void exploreCell(meshContainer_t *meshContainer, cell_t *currentCell, int nStepFromStart, cell_t *startCell, cell_t *endCell) {
   // static unsigned int bestLen = 999999;
    offsets_t moves[26];
#ifdef DEBUG
    fprintf(stderr, "%d %d %d %d\n", currentCell->i, currentCell->j, currentCell->k, currentCell->memberCount);
#endif
    printf("%d %d %d %d\n", currentCell->i, currentCell->j, currentCell->k, currentCell->bwfs);
    if (areSameCells(currentCell, endCell)) {
        printf("Found destination at depth %d\n", nStepFromStart);
        bestLen = nStepFromStart;
    }
    if (nStepFromStart > 0 && areSameCells(currentCell, startCell) ) {
        printf("BACK TO SQ ONE\n");
        return;
    }
    if (nStepFromStart > 0 && currentCell->memberCount > 0)
        return; // Cell is blocked

    if (nStepFromStart > currentCell->bwfs)
        return;
    currentCell->bwfs = nStepFromStart;
    if (c_dist(currentCell,  endCell) + currentCell->bwfs > bestLen) {
        fprintf(stderr, "%d %d %d to far (%f)\n", \
        currentCell->i, currentCell->j, currentCell->k,\
        c_dist(currentCell,  endCell) + currentCell->bwfs);
        return;
    }
    
    short int neighbourCount = sortNeighboursByMeshDistance(currentCell, endCell, meshContainer->mesh, moves);
    for (int iCell = 0 ; iCell < neighbourCount ; iCell++) {
        int i = moves[iCell].i;
        int j = moves[iCell].j;        
        int k = moves[iCell].k;
        exploreCell(meshContainer, &meshContainer->mesh->grid[i][j][k],\
                    nStepFromStart + 1, startCell, endCell);
    }

}

float c_dist(cell_t *a, cell_t *b) {
    return sqrt(  ( a->i - b->i ) * ( a->i - b->i ) \
                + ( a->j - b->j ) * ( a->j - b->j ) \
                + ( a->k - b->k ) * ( a->k - b->k ) \
            );
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
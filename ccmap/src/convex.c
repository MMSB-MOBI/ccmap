#include "convex.h"
// Overflow of grid is possible ^^
// We need max radius + one cell for path 
int buildSphere(atom_t *atom, cell_t *optCell, meshContainer_t *meshContainer) {
    // Compute radius of current atom in cell units
    int rad_cu = atom->_radiusASA / meshContainer->step + 0.5; // should be rounded up
    // Get corresponding cell and compute diameter boundaries index along the 3 axis
    int rac_cu_pow = rad_cu * rad_cu;
    cell_t *cCell = optCell == NULL ?\
        getCellFromAtom(meshContainer, atom):\
        optCell;
//#ifdef DEBUG
    char bufferLog[81];
    stringifyAtom(atom, bufferLog);
    fprintf(stderr, "Building voxels sphere for:%s\n", bufferLog);
//#endif
    // filled up cells
    int nvx = 0;
    cell_t *currCell = NULL;
    for (int i = cCell->i - rad_cu ; i <= cCell->i + rad_cu ; i++) {
        for (int j = cCell->j - rad_cu ; j <= cCell->j + rad_cu; j++) {
            for (int k = cCell->k - rad_cu ; k <= cCell->k + rad_cu; k++) {   
                currCell = &meshContainer->mesh->grid[i][j][k];
                int norm = i * i + j * j + k * k;
                if (norm > rac_cu_pow)
                    continue;
                nvx++;
                if (norm == rac_cu_pow && !currCell-> isInterior)
                    currCell->isSurface = true;
                if(norm < rac_cu_pow) {
                    if(currCell->isSurface)
                        nvx--; // voxel was already accounted for
                    currCell->isInterior = true;
                    currCell->isSurface  = false;
                }
            }
        }
    }
    return nvx; // overestimation as surface voxel downgrade as volume are counted
}

// Custom predicate to allow/reject cell exploration in pathfinder
bool surfaceExplorerPredicate(cell_t *cell){
    return cell->memberCount == 0 || cell->isSurface;
}

void buildSurfaces(meshContainer_t *meshContainer) {
    int totalVoxel = 0;
    cell_t *currCell = NULL;
    // iterate over each non empty cell
    for(int i_cell = 0 ; i_cell < meshContainer->nFilled ; i_cell++) {
        currCell = meshContainer->filledCells[i_cell];
        if (currCell->memberCount > 1) {
            fprintf(stderr, \
            "Fatal:buildSurface cell %d %d %d does not hold one single atom (%d)\n",\
            currCell->i, currCell->j, currCell->k, currCell->memberCount\
            );
            return -1;
        }
        // members field is a straight pointer --> Head of a chained list        
        totalVoxel += buildSphere(currCell->members, currCell,\
                            meshContainer);
    }
    // assert only one atom at most
    // build sphere
    // summing filled cells
    meshContainer->nVoxels = totalVoxel;
}

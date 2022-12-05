#include "convex.h"

// Custom predicate to allow/reject cell exploration in pathfinder
bool surfaceExplorerPredicate(cell_t *cell){
    return (!cell->isInterior) && (cell->memberCount == 0 || cell->isSurface);
}

bool buildSurfaces(meshContainer_t *meshContainer) {
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
            return false;
        }
        // members field is a straight pointer --> Head of a chained list        
        totalVoxel += buildSphere(currCell->members, currCell,\
                            meshContainer);
    }
    // assert only one atom at most
    // build sphere
    // summing filled cells
    meshContainer->nVoxels = totalVoxel;
    return true;
}

setCells_t *getSurfaceCells(atom_t *atom, meshContainer_t *meshContainer) {
    int rad_cu = atom->_radiusASA / meshContainer->step + 0.5;
    
    cell_t *cCell = getCellFromAtom(meshContainer, atom);
    cell_t *currCell = NULL;
    setCells_t *setCells = malloc(sizeof(setCells_t));
    setCells->size = 0;
    for (int i = cCell->i - rad_cu ; i <= cCell->i + rad_cu ; i++) {
        for (int j = cCell->j - rad_cu ; j <= cCell->j + rad_cu; j++) {
            for (int k = cCell->k - rad_cu ; k <= cCell->k + rad_cu; k++) { 
                currCell = &meshContainer->mesh->grid[i][j][k];
                if(currCell->isSurface && ! currCell->isInterior)
                    setCells->size++;
    }}}
    setCells->cells = malloc( setCells->size * sizeof(cell_t*) );
    int i_cell = 0;
    for (int i = cCell->i - rad_cu ; i <= cCell->i + rad_cu ; i++) {
        for (int j = cCell->j - rad_cu ; j <= cCell->j + rad_cu; j++) {
            for (int k = cCell->k - rad_cu ; k <= cCell->k + rad_cu; k++) { 
                currCell = &meshContainer->mesh->grid[i][j][k];
                if(currCell->isSurface && ! currCell->isInterior) {
                    setCells->cells[i_cell] = currCell;
                    i_cell++; 
                }
    }}}
    return setCells;
}

int buildSphere(atom_t *atom, cell_t *optCell, meshContainer_t *meshContainer) {
    // Compute radius of current atom in cell units
    int rad_cu = atom->_radiusASA / meshContainer->step + 0.5; // should be rounded up
    // Get corresponding cell and compute diameter boundaries index along the 3 axis
    int rad_cu_pow = rad_cu * rad_cu;
    cell_t *cCell = optCell == NULL ?\
        getCellFromAtom(meshContainer, atom):\
        optCell;
#ifdef DEBUG
    char bufferLog[81];
    stringifyAtom(atom, bufferLog);
    fprintf(stderr, "Building voxels sphere(rad=%d around %d %d %d)\n\t=>aka:%s\n",\
        rad_cu, cCell->i, cCell->j, cCell->k, bufferLog);
    
#endif
    // filled up cells
    int nvx = 0;
    cell_t *currCell = NULL;
    for (int i = cCell->i - rad_cu ; i <= cCell->i + rad_cu ; i++) {
        for (int j = cCell->j - rad_cu ; j <= cCell->j + rad_cu; j++) {
            for (int k = cCell->k - rad_cu ; k <= cCell->k + rad_cu; k++) {   
                currCell = &meshContainer->mesh->grid[i][j][k];
                nvx += voxelEvaluate(currCell, cCell, (double)rad_cu_pow);
    }}}
    return nvx;// voxel actually constructed
}

// We stick to cell space
// compute distance between 8 corners of curr_cell and center of center cell
int voxelEvaluate(cell_t *currCell, cell_t *centerCell, double norm) {
    double i0 = ((double)centerCell->i) + 0.5;
    double j0 = ((double)centerCell->j) + 0.5;
    double k0 = ((double)centerCell->k) + 0.5;
    int insideCorners = 0;
    double currNorm;
    // Loop over current cell corners
    for (int iOff = 0 ; iOff <=1 ; iOff++) {
        for (int jOff = 0 ; jOff <=1 ; jOff++) {
            for (int kOff = 0 ; kOff <=1 ; kOff++) {
                double iC = (double) (currCell->i + iOff);
                double jC = (double) (currCell->j + jOff);
                double kC = (double) (currCell->k + kOff);

                currNorm = (i0 - iC) * (i0 - iC)\
                         + (j0 - jC) * (j0 - jC)\
                         + (k0 - kC) * (k0 - kC);
                if (currNorm <= norm)
                    insideCorners++;
            }
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Voxel (%d %d %d) buried corners count %d\n",\
            currCell->i, currCell->j, currCell->k, insideCorners   
        );
#endif
   
    if (insideCorners == 0) // Outside current sphere
        return 0;
    if (insideCorners == 8) { //it is a buried voxel
        if (currCell->isInterior) // already registred as filed
            return 0;
        currCell->isInterior = true;
        if (currCell->isSurface) // was a surface of another sphere, dont increment total volume voxel count
            return 0;
        return 1;// Declare new volume voxel
    }
    // it is a surface voxel
    if (currCell->isSurface) // was a surface of another sphere, dont increment total volume voxel count
        return 0;
    currCell->isSurface = true;
    if (currCell->isInterior) // was a buried voxel of another sphere, dont increment total volume voxel count
        return 0;
        
    return 1; // Declare new volume voxel
}

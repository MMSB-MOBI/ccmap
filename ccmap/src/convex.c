#include "convex.h"

// Custom predicate to allow/reject cell exploration in pathfinder
bool surfaceExplorerPredicate(cell_t *cell){
    return (!cell->isInterior) && (cell->memberCount == 0 || cell->isSurface);
}

bool buildSurfaces(meshContainer_t *meshContainer, bool force) {
    int totalVoxel = 0;
    cell_t *currCell = NULL;
    char buffer[81];
    struct atom *buffAtom = NULL;
    // iterate over each non empty cell
    for(int i_cell = 0 ; i_cell < meshContainer->nFilled ; i_cell++) {
        currCell = meshContainer->filledCells[i_cell];
        if (currCell->memberCount > 1) {
            
            buffAtom = currCell->members;
            char *errType = force ? "Warning" : "Fatal";
            fprintf(stderr, \
            "%s:buildSurface cell %d %d %d does not hold one single atom (%d)\n",\
            errType, currCell->i, currCell->j, currCell->k, currCell->memberCount\
            );
           
            for(int i = 0 ; i < currCell->memberCount ; i++) {
                stringifyAtom(buffAtom, buffer);
                fprintf(stderr, "[member %d] %s\n", i + 1, buffer);
                buffAtom = buffAtom->nextCellAtom;
            }
            if(!force)
                return false;  
        }
        // members field is a straight pointer --> Head of a chained list        
        totalVoxel += buildSphere(currCell->members, currCell,\
                            meshContainer);
    }

    meshContainer->nVoxels = totalVoxel;
    return true;
}


/*
    Get all voxels which are "surface" of prvided atom, as a set.
*/
setCells_t *getSurfaceCells(atom_t *atom, meshContainer_t *meshContainer) {

    int rad_cu = atom->_radiusASA / meshContainer->step + 0.5;
    

    /*
    bool dumpMe = false; 
    if (strcmp( atom->resName, "LI1") == 0 && strcmp( atom->name, " CA ")==0 )
        dumpMe = true;
    */


  /*  cell_t *cCell = getCellFromAtom(meshContainer, atom);
    cell_t *currCell = NULL;
    setCells_t *setCells = malloc(sizeof(setCells_t));
    setCells->size = 0;
    bool isValidCell = false;
    static offsets[] = { -1 , 1 };
    double cx, cy, cz;
    */
    setCells_t *setCells = malloc(sizeof(setCells_t));
    setCells->size = 0;
    cell_t *cCell = getCellFromAtom(meshContainer, atom);
    cell_t *currCell = NULL;
    int interiorCorners;
    for (int i = cCell->i - rad_cu ; i <= cCell->i + rad_cu ; i++) {
        for (int j = cCell->j - rad_cu ; j <= cCell->j + rad_cu; j++) {
            for (int k = cCell->k - rad_cu ; k <= cCell->k + rad_cu; k++) { 
                currCell = &meshContainer->mesh->grid[i][j][k];
                if(!currCell->isSurface || currCell->isInterior) 
                    continue;
                
                // now check that this surface voxel belongs to selection atom                
                interiorCorners = countCorners(meshContainer, atom, currCell);
                if(interiorCorners == 0){
                  //  if(dumpMe)  
                   //     printf("!!! Voxel %d %d %d was a false positive!!\n", currCell->i, currCell->j, currCell->k);
                    continue;
                }
                setCells->size++;
    }}}
    setCells->cells = malloc( setCells->size * sizeof(cell_t*) );
    int i_cell = 0;
    for (int i = cCell->i - rad_cu ; i <= cCell->i + rad_cu ; i++) {
        for (int j = cCell->j - rad_cu ; j <= cCell->j + rad_cu; j++) {
            for (int k = cCell->k - rad_cu ; k <= cCell->k + rad_cu; k++) { 
                currCell = &meshContainer->mesh->grid[i][j][k];

                if(!currCell->isSurface || currCell->isInterior) 
                    continue;
                
                // now check that this surface voxel belongs to selection atom                
                interiorCorners = countCorners(meshContainer, atom, currCell);
                if(interiorCorners == 0) 
                    continue;
                setCells->cells[i_cell] = currCell;
                i_cell++; 

    }}}
    return setCells;
}

/*
    Count the number of buried corners with of cell inside atom exclusion sphere
*/
int countCorners(meshContainer_t *meshContainer, atom_t *atom, cell_t *cell) {

    double radius = atom->_radiusASA;
    double step   = meshContainer->step / 2;
    static int offsets[] = { -1 , 1 };
    int insideCorners = 0;
    double _x, _y, _z;
    double d;
    // Explicit distance computation of all cube vertex eucldiean dist w/ center atom

    double cx, cy, cz;
    meshToCartesian(meshContainer, cell->i, cell->j, cell->k, &cx, &cy, &cz);


   /*
   bool dumpMe = false; 
    if (strcmp( atom->resName, "LI1") == 0 && strcmp( atom->name, " CA ")==0 )
        dumpMe = true;
    */
    // Loop over current cell corners
    for (int i = 0 ; i <=1 ; i++) {
        for (int j = 0 ; j <=1 ; j++) {
            for (int k = 0 ; k <=1 ; k++) {
                _x = cx + offsets[i] * step;
                _y = cy + offsets[j] * step;
                _z = cz + offsets[k] * step;
                d = euclideanDistance3(atom->x, atom->y, atom->z,\
                                       _x, _y, _z);
                /*
                if(dumpMe)
                    printf("curr_center %g %g %g --> offsets %d %d %d --> corners %g %g %g VS atom %g %g %g => d:%g\n",\
                        cx, cy, cz, i, j, k, _x, _y, _z, atom->x, atom->y, atom->z, d);
                */
                if(d <= radius)
                    insideCorners++;
            }
        }
    }
    return insideCorners;
}

int buildSphere(atom_t *atom, cell_t *optCell, meshContainer_t *meshContainer) {
    // Compute radius of current atom in cell units
    int rad_cu = atom->_radiusASA / meshContainer->step + 0.5; // should be rounded up
    // Get corresponding cell and compute diameter boundaries index along the 3 axis
   // int rad_cu_pow = rad_cu * rad_cu;
    cell_t *cCell = optCell == NULL ?\
        getCellFromAtom(meshContainer, atom):\
        optCell;
#ifdef DEBUG
    char bufferLog[81];
    stringifyAtom(atom, bufferLog);
    fprintf(stderr, "Building voxels sphere(rad=%d around %d %d %d)\n\t=>aka:%s\n",\
        rad_cu, cCell->i, cCell->j, cCell->k, bufferLog);
    
#endif
   /*
    bool dumpMe = false; 
    if (strcmp( atom->resName, "LI1") == 0 && strcmp( atom->name, " CA ")==0 )
        dumpMe = true;
    
    if(dumpMe) {
        char bufferLog[81];
        stringifyAtom(atom, bufferLog);
        fprintf(stderr, "Building voxels sphere(rad=%d around %d %d %d)\n\t=>aka:%s\n",\
            rad_cu, cCell->i, cCell->j, cCell->k, bufferLog);
    }
    */
    // filled up cells
    /*
    double _x, _y, _z;
    double _cx, _cy, _cz;
    double _d, _d2;
    */
    int nvx = 0;
    cell_t *currCell = NULL;
    for (int i = cCell->i - rad_cu ; i <= cCell->i + rad_cu ; i++) {
        for (int j = cCell->j - rad_cu ; j <= cCell->j + rad_cu; j++) {
            for (int k = cCell->k - rad_cu ; k <= cCell->k + rad_cu; k++) {   
                currCell = &meshContainer->mesh->grid[i][j][k];
                //nvx += voxelEvaluate(currCell, cCell, (double)rad_cu_pow);
                nvx += voxelEvaluateCartesian(meshContainer, currCell, atom);
                /*
                if(currCell->isSurface dumpMe) {
                    meshToCartesian(meshContainer, i, j, k, &_x, &_y, &_z);
                    meshToCartesian(meshContainer, cCell->i, cCell->j, cCell->k, &_cx, &_cy, &_cz);
                    _d = euclideanDistance3(_x, _y, _z, _cx, _cy, _cz);
                    printf("RAD_CU  %d = %g / %g +0.5\n", rad_cu, atom->_radiusASA, meshContainer->step);
                    printf("Current Surf Cell [%d %d %d // %g %g %g ] cart distance from center cell [%d %d %d // %g %g %g] : %g\n",\
                        i, j, k, _x, _y,_z, cCell->i,cCell->j, cCell->k, _cx, _cy, _cz, _d);
                    _d2 = euclideanDistance3(_x, _y, _z, atom->x, atom->y, atom->z);
                    printf("Current Surf Cell [%d %d %d // %g %g %g ]cart distance from center atom [%d %d %d // %g %g %g] : %g\n",\
                        i, j, k, _x, _y,_z, cCell->i,cCell->j, cCell->k, atom->x, atom->y, atom->z, _d2);
                    //printf("delta: %g\n", fabs(_d - _d2));
                }
                */

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

// We define surface and buried voxel based on cartesian metric from actual atom center
int voxelEvaluateCartesian(meshContainer_t *meshContainer, cell_t *currCell, atom_t *centerAtom) {
    
    int insideCorners = countCorners(meshContainer, centerAtom, currCell);

    /*
    bool dumpMe = false; 
    if (strcmp( centerAtom->resName, "LI1") == 0 && strcmp( centerAtom->name, " CA ")==0 )
        dumpMe = true;
*/
#ifdef DEBUG
    fprintf(stderr, "Voxel (%d %d %d) buried corners count %d\n",\
            currCell->i, currCell->j, currCell->k, insideCorners   
        );
#endif
  /*  if(dumpMe)
        printf("Voxel (%d %d %d) buried corners count %d\n\n",\
            currCell->i, currCell->j, currCell->k, insideCorners   
        );
*/
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

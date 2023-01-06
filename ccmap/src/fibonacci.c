#include "fibonacci.h"


// Compute Sphere Total and buried ASA
// We compute by decreasing total
void computeFiboSphereASA(fibo_grid_t *iFiboGrid, float *totalSurface, float *buriedSurface) {

    *(totalSurface)   = (float)iFiboGrid->radius * (float)iFiboGrid->radius * (float)4.00 * M_PI;


  //  fprintf(stderr, "%d\n",iFiboGrid->n_spots);
    float spotSurface = *(totalSurface) / (float)iFiboGrid->n_spots;
    *(buriedSurface)  = 0;
    #ifdef DEBUG
        int nbBuried = 0;

        fprintf(stderr,"Compute FiboSphere (r:%f, pi:%f) Total/Spot %f %f\n",\
        iFiboGrid->radius, M_PI, *totalSurface, spotSurface);
    #endif
    for (int k = 0; k < iFiboGrid->n_spots ; k++) {
        #ifdef DEBUG
            fprintf(stderr,"=FIBO= %d * %f\n", iFiboGrid->spots[k].buried, spotSurface);
            nbBuried += iFiboGrid->spots[k].buried;
        #endif
        *(buriedSurface) += (float)iFiboGrid->spots[k].buried * spotSurface;
    }
    #ifdef DEBUG
        fprintf(stderr, "Computing FiboSphereASA rad = %f total = %f  buried %f [unit:%f, nb_buried:%d, total:%d]\n",\
                iFiboGrid->radius, *totalSurface, *buriedSurface, spotSurface, nbBuried, FIBO_K);
    #endif
    
    // Rounding up and summing can exceed exact total
    *(buriedSurface) = *buriedSurface > *totalSurface ? *totalSurface : *buriedSurface;
}

// Computes if i grid spots penetrates j sphere// we perform a similar operation on  j grid spots
void FiboSpherePairProcess(fibo_grid_t *iFiboGrid, fibo_grid_t *jFiboGrid) {
    #ifdef DEBUG
        int nb_buried_iSpots = 0;
        int nb_buried_jSpots = 0;
    #endif
    assert(iFiboGrid->n_spots == jFiboGrid->n_spots);
    for (int k = 0; k < iFiboGrid->n_spots ; k++) {
        if( !iFiboGrid->spots[k].buried ){ //continue;
            iFiboGrid->spots[k].buried = \
                euclideanDistance3(\
                    iFiboGrid->spots[k].x, iFiboGrid->spots[k].y, iFiboGrid->spots[k].z,\
                    jFiboGrid->center.x  , jFiboGrid->center.y  , jFiboGrid->center.z)\
                < jFiboGrid->radius;
        }
        if( !jFiboGrid->spots[k].buried ){ //continue;
            jFiboGrid->spots[k].buried = \
                euclideanDistance3(\
                    jFiboGrid->spots[k].x, jFiboGrid->spots[k].y, jFiboGrid->spots[k].z,\
                    iFiboGrid->center.x  , iFiboGrid->center.y  , iFiboGrid->center.z)\
                < iFiboGrid->radius;
        }
        #ifdef DEBUG
            nb_buried_iSpots += iFiboGrid->spots[k].buried;
            nb_buried_jSpots += jFiboGrid->spots[k].buried;
        #endif
    }
    #ifdef DEBUG
        fprintf(stderr, "CURRENT TOTAL_BURIED_SPOT i/j = %d\n", nb_buried_iSpots, nb_buried_jSpots);
    #endif
}


fibo_grid_t *computeFiboGrid(float x, float y, float z, float radius, bool hres) {
    int FIBO_K = hres ? FIBO_K20 : FIBO_K14;
    int FIBO_K_min1 = hres ? FIBO_K19 : FIBO_K13;
    
    #ifdef DEBUG
    fprintf(stderr, "Starting atomicFiboGrid [%g, %g, %g] r=%f\n", x, y, z, radius);
    #endif
    //radius = 50;
  
    fibo_grid_t *fibo_grid = createFiboGrid(FIBO_K);
    
    int phi_prime = 0;
    float theta, phi;
    fibo_grid->center.x = x;
    fibo_grid->center.y = y;
    fibo_grid->center.z = z;
    fibo_grid->radius   = radius;

    for (int i_spot = 0 ; i_spot < FIBO_K ; i_spot++) {
        phi_prime = phi_prime + FIBO_K_min1 <= FIBO_K ? \
                    phi_prime + FIBO_K_min1 :\
                    phi_prime + FIBO_K_min1 - FIBO_K;

		theta = acos(1.0-2.0 * i_spot / FIBO_K);
		phi   = 2.0 * M_PI * phi_prime / FIBO_K;
		fibo_grid->spots[i_spot].x = x + radius * sin(theta)*cos(phi);
        fibo_grid->spots[i_spot].y = y + radius * sin(theta)*sin(phi);
        fibo_grid->spots[i_spot].z = z + radius * cos(theta);
        fibo_grid->spots[i_spot].buried = false;
    }

    return fibo_grid;
}

void updateFiboGrid(fibo_grid_t *fibo_grid, float dx, float dy, float dz) {
    fibo_grid->center.x += dx;
    fibo_grid->center.y += dy;
    fibo_grid->center.z += dz;
    for (int i_spot = 0 ; i_spot < fibo_grid->n_spots ; i_spot++) {
        fibo_grid->spots[i_spot].x += dx;
        fibo_grid->spots[i_spot].y += dy;
        fibo_grid->spots[i_spot].z += dz;
        fibo_grid->spots[i_spot].buried = false;
    }
}
/*
void updateFiboGrid(fibo_grid_t *fibo_grid, float x, float y, float z, bool hres) {
    int FIBO_K = hres ? FIBO_K20 : FIBO_K14;
    int FIBO_K_min1 = hres ? FIBO_K19 : FIBO_K13;
    #ifdef DEBUG
    fprintf(stderr, "Starting update atomicFiboGrid [%g, %g, %g] r=%f\n", x, y, z, radius);
    #endif
  
    
    int phi_prime = 0;
    float theta, phi;
    fibo_grid->center.x = x;
    fibo_grid->center.y = y;
    fibo_grid->center.z = z;
    // fibo_grid->radius; Already set
    float radius = fibo_grid->radius;

    for (int i_spot = 0 ; i_spot < FIBO_K ; i_spot++) {
        phi_prime = phi_prime + FIBO_K_min1 <= FIBO_K ? \
                    phi_prime + FIBO_K_min1 :\
                    phi_prime + FIBO_K_min1 - FIBO_K;

		theta = acos(1.0-2.0 * i_spot / FIBO_K);
		phi   = 2.0 * M_PI * phi_prime / FIBO_K;
		fibo_grid->spots[i_spot].x = x + radius * sin(theta)*cos(phi);
        fibo_grid->spots[i_spot].y = y + radius * sin(theta)*sin(phi);
        fibo_grid->spots[i_spot].z = z + radius * cos(theta);
        fibo_grid->spots[i_spot].buried = false;
    }
    
}
*/



fibo_grid_t *createFiboGrid(int n) {
    fibo_grid_t *fibo_grid = malloc( sizeof(fibo_grid_t) );
    fibo_grid->spots = malloc( sizeof(spot_t) * n );
    fibo_grid->n_spots = n;

    return fibo_grid;
}

fibo_grid_t *destroyFiboGrid(fibo_grid_t *fibo_grid){
    free(fibo_grid->spots);
    free(fibo_grid);

    return fibo_grid;
}

string_t *jsonifyFiboGrid(fibo_grid_t *fibo_grid){

    #ifdef DEBUG
    fprintf(stderr, "Starting jsonifyFiboGrid\n");
    #endif
    string_t *jsonString = createString();

    if (fibo_grid == NULL) {
        jsonString->append(jsonString, "{}");
        return jsonString;
    }
    char buffer[1024];
    sprintf(buffer, "{\"center\": [%.3g, %.3g, %.3g], \"grid\":[", \
        fibo_grid->center.x, fibo_grid->center.y, fibo_grid->center.z);
    jsonString->append(jsonString, buffer);
    for (int i = 0 ; i < fibo_grid->n_spots ; i++) {
        sprintf(buffer,"[%.3f, %.3f, %.3f]", \
            fibo_grid->spots[i].x, \
            fibo_grid->spots[i].y, \
            fibo_grid->spots[i].z);

        jsonString->append(jsonString, buffer);
        if(i < (fibo_grid->n_spots - 1))
            jsonString->append(jsonString, ",\n");

    }
    jsonString->append(jsonString, "]}");
   
    return jsonString;
}

void printFiboGrid(fibo_grid_t *fibo_grid){
    string_t *stringifyFigoGrid = jsonifyFiboGrid(fibo_grid);
    char *s = stringifyFigoGrid->toChar(stringifyFigoGrid);
    printf( "%s", s );
    free(s);
    destroyString(stringifyFigoGrid);
}

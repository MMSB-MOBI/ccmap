#include "fibonacci.h"

fibo_grid_t *computeFiboGrid(float x, float y, float z, float radius) {
    #ifdef DEBUG
    fprintf(stderr, "Starting atomicFiboGrid\n");
    #endif
    //radius = 50;
    fibo_grid_t *fibo_grid = createFiboGrid(FIBO_K14);
    
    int phi_prime = 0;
    float theta, phi;
    fibo_grid->center.x = x;
    fibo_grid->center.y = y;
    fibo_grid->center.z = z;

    for (int i_spot = 0 ; i_spot < FIBO_K14 ; i_spot++) {
        phi_prime = phi_prime + FIBO_K13 <= FIBO_K14 ? \
                    phi_prime + FIBO_K13 :\
                    phi_prime + FIBO_K13 - FIBO_K14;

		theta = acos(1.0-2.0 * i_spot / FIBO_K14);
		phi   = 2.0 * M_PI * phi_prime / FIBO_K14;
		fibo_grid->spots[i_spot].x = x + radius * sin(theta)*cos(phi);
        fibo_grid->spots[i_spot].y = y + radius * sin(theta)*sin(phi);
        fibo_grid->spots[i_spot].z = z + radius * cos(theta);
        fibo_grid->spots[i_spot].buried = false;
    }

    return fibo_grid;
}

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
        if(i < fibo_grid->n_spots - 1)
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

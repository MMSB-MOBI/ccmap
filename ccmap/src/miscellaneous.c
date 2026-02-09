#include "miscellaneous.h"

double euclideanDistance3(float x0, float y0, float z0, float x1, float y1, float z1) {
    return sqrt( (x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1) );
}

#ifdef DEBUG
void printOnContextStderr(char *debugCharContent) {
        #ifdef AS_PYTHON_EXTENSION
            PySys_WriteStderr("%s", debugCharContent);
        #else
            fprintf(stderr, "%s", debugCharContent);
        #endif
}
#endif

void strip(char *dest, char *src) { // Remove space from src into dest
    int i,j;
    i = 0;
    j = 0;
    while(src[j] != '\0') {
        if (src[j] != '\t' && src[j] != ' ') {
            dest[i] = src[j];
            i++;
        }
        j++;
    }
    dest[i] = '\0';

}


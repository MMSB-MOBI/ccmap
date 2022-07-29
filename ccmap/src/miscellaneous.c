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

#include "parameters.h"

void stringToThreeFloats(char *input, float (*vector)[3]) {
    char *err, *p = input;
    if (input == NULL) return;
    float val;
    int i = 0;
    while (*p) {
        val = strtod(p, &err);
        if (p == err) p++;
        else if ((err == NULL) || (*err == 0)) {
            (*vector)[i] = val;
            i++;
        //    printf("Value: %f\n", val);
            break;
        }
        else {
      //      printf("errValue: %f\n", val);
            p = err + 1;
            (*vector)[i] = val;
            i++;
        }
    }
    //printf("->%g %g %g\n", (*vector)[0], (*vector)[1], (*vector)[2]);

}


void parseTransform(char *eulerString, char *translateString, float (*eulers)[3], float (*translate)[3]){
    if (translateString != NULL){
        stringToThreeFloats(translateString, translate);
        //printf ("Cartesian translation vector %g %g %g\n", (*translate)[0], (*translate)[1], (*translate)[2]);
    } else {
       // translate = NULL;
    }
    if (eulerString != NULL){
        stringToThreeFloats(eulerString, eulers);
        //printf ("Euler's angles values %g %g %g\n", (*eulers)[0], (*eulers)[1], (*eulers)[2]);
    } else {
        //eulers = NULL;
    }

}

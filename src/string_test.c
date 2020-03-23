#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "default.h"


int main () {
    char hello[] = "Hello"; // Will include '\0'
    char morning[] = " DELU";

    string_t *s = createString();
    s->dump(s);
    s->append(s, hello);
    s->dump(s);
    s->append(s, morning);

    printf("[%d,%lu]\"%s\"\n", s->length, s->capacity, s->value);
    s->dump(s);

    string_t *s2 = s->copy(s);
    s2->dump(s2);


    destroyString(s);
    destroyString(s2);
    return 1;
}

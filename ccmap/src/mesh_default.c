#include "mesh_default.h"
/*
strcpy Copies the C string pointed by source into the array
pointed by destination, including the terminating null character
(and stopping at that point).
strlen(const char *str) computes the length
of the string str up to, but not including the terminating null character.
*/
// Both chains are supposed to have a teminating null char Alternatively dest can be a NULL string pointer
// Concatenate two strings in the 1st one, return the length of the string not including the '\0' termination character
int concatenate(char **dest, char *src) {
    int oldStringLength = *dest != NULL ? strlen(*dest) : 0;
    // oldStrinLength is the index at which we will start to append/copy the src srting into dest string
    int new_buf_size = oldStringLength + strlen(src) + 1; // Both string length + slot for '\0'

    *dest = realloc( *dest, new_buf_size * sizeof(char) );
    strcpy(&((*dest)[oldStringLength]), src);
    if ((*dest)[new_buf_size - 1] != '\0') {
        printf("String copy buffer termination error");
    }
    // This line should not be necessary
    (*dest)[new_buf_size - 1] = '\0';

    return new_buf_size - 1;
}

// Pop the last element of a string if it matches provided char c
// The length of the string is effectively reduced by one, this is the returned value
int popChar(char **dest, char c) {
    int i = 0;
    while( (*dest)[i] != '\0' ) {
        i++;
    }
    if ( (*dest)[i - 1] == c ) {
        (*dest)[i - 1] = '\0';
         *dest = realloc( *dest, (i) * sizeof(char) );
    }
    return i;
}


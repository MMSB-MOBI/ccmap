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

void dumpString (string_t* this) {
    printf("%d\n", this->length);
    for (int i = 0; i < this->length + 1; i++) {
        printf("%d:\"%c\"\n", i, this->value[i]);
    }
}
bool appendString(string_t *this, char *word) {
    size_t wLen = strlen(word);
    if (this->length + wLen >= this->capacity) { //>= to rpeserve 1 slot for '\0'
     //   fprintf(stderr, "Extending string from %lu to %lu\n", this->capacity, this->capacity + STRING_CHUNCK);
        this->capacity += STRING_CHUNCK;
        char *_value = realloc(this->value, this->capacity);
        assert(_value != NULL);
        this->value = _value;
        this->append(this, word);
        return true;
    }
    strcpy( &(this->value[this->length]), word );
    this->length += wLen; 
    this->value[this->length] = '\0';
    return true;
}

string_t *copyString(string_t *this) {
    string_t *other = createString();
    other->append(other, this->value);

    return other;
}

string_t *createString() {
    string_t *this = malloc(sizeof(string_t));
    this->capacity = STRING_CHUNCK;
    this->value = malloc(this->capacity * sizeof(char));
    this->value[0] = '\0';
    this->length = 0;
    this->append = &appendString;
    this->copy = &copyString;
    this->dump = &dumpString;
    this->toChar = &toCharString;
    return this;
} 

char *toCharString(string_t *this) {
    char *asChar = malloc((this->length + 1) * sizeof(char));
    strcpy(asChar, this->value);
    return asChar;
}
string_t *destroyString(string_t *this) {
    free(this->value);
    free(this);
    return this;
}

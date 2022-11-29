#include "my_string.h"

// Return datastructure of a split method

void dumpString (string_t* this) {
    //printf("%d\n", this->length - 1);
    /*for (int i = 0; i < this->length + 1; i++) {
        printf("%d:\"%c\"\n", i, this->value[i]);
    }*/
    printf("%d:\"%s\"\n", this->length, this->value);
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

// Split into strings management
stringList_t *splitAndCreateStringList(char *input, char sep) {
    stringList_t *stringList = malloc( sizeof(stringList_t) );
    string_t *string;
    char buffer[81];
    int i_char = 0;
    int i = 0;
    int i_field = 1;
    while (input[i_char] != '\0') {
        if (input[i_char] == sep) 
            i_field++;
        i_char++;
    }
    /*
    fprintf(stderr, "%s==>%d\n", input, i_field);
    */
    stringList->elem = malloc( i_field * sizeof(string_t*) );
    stringList->len = i_field;
    i_char  = 0;
    i_field = 0;
    while (input[i_char] != '\0') {
        if (input[i_char] == sep){
            buffer[i_field] = '\0';
            stringList->elem[i] = createString();
            string = stringList->elem[i];         
            string->append(string, buffer);

        /*
            fprintf(stderr, "Putting \"%s\" at %d\n", buffer, i);
            string->dump(string);
        */
            i_char++;
            i_field = 0;
            i++;
            continue;
        }
        
        

        buffer[i_field] = input[i_char];
        i_field++;
        i_char++;
    }
    stringList->elem[i] = createString();
    string = stringList->elem[i];
    buffer[i_field] = '\0';
    string->append(string, buffer);
    /*
    fprintf(stderr, "Putting \"%s\" at %d\n", buffer, i);
    string->dump(string);
    */
    return stringList;
}

stringList_t *destroyStringList(stringList_t *stringList) {
    for (int i = 0 ; i < stringList->len ; i++)
        stringList->elem[i] = destroyString(stringList->elem[i]);
    free(stringList->elem);
    free(stringList);
    return NULL;
}

bool compareStringToChar(string_t *self, char *other, bool bStrip) {
    char buff_i[81]; 
    char buff_j[81];
    char *stringBuf = toCharString(self);
    if (!bStrip)
        return strcmp(other, stringBuf) == 0;   

    strip(buff_i, other);
    strip(buff_j, stringBuf);
    free(stringBuf);
    return strcmp(buff_i, buff_j) == 0;    
}
#include "atom_mapper.h"

void create_buffers(char ***names_buffer, char ***residue_names_buffer, float **radii_buffer){
    *names_buffer         = malloc( MAX_ATOM_PAYLOAD * sizeof(char**) );
    *residue_names_buffer = malloc( MAX_ATOM_PAYLOAD * sizeof(char**) );
    *radii_buffer         = malloc ( MAX_ATOM_PAYLOAD * sizeof(float) );

}
/* not sure /bt types here
void destroy_buffer(char ***names_buffer, char ***residue_names_buffer, float **radii_buffer){
    // TO DO
}
*/

atom_map_t *readAtomMapperFromFile(char *filePath, float probeRadius) {
    fprintf(stderr, "readAtomMapperFromFile\n");
    fprintf(stderr, "==>%s\n", filePath);
    FILE *fp;
    atom_map_t *atom_map;
    int   cur_count;
    float cur_radius;
    char  cur_name[81];
    char  cur_resname[81];
    char **atom_name_buffer = NULL;
    char **res_name_buffer  = NULL;
    float *radii_buffer     = NULL;
    create_buffers(&atom_name_buffer, &res_name_buffer, &radii_buffer);

    fprintf(stderr, "toto\n");
    cur_name[0] = '\0';
    if ( !(fp = fopen(filePath, "r")) ) {
        fprintf(stderr, "File open error at %s\n", filePath);
        return NULL;
    }
    fprintf(stderr, "toto\n");
    while (fscanf(fp, "%s\t%s\t%f\n", cur_resname, cur_name, &cur_radius) != EOF) {
        printf(">%s\t%s\t%f\n", cur_resname, cur_name, cur_radius);
    }

    fclose(fp);
    return atom_map;
}

atom_map_t *destroyAtomMapper(atom_map_t *map){
    return NULL;
}

/*
atom_map_t *readAtomMapperFromFile(char *filePath, float probeRadius);
void readMapperInputs(FILE *fp, char **atom_names, char **residue_names, float **radii, float probeRadius);
atom_map_t *createAtomMapper(char **atom_names, char **residue_names, float **radii);
atom_map_t *destroyAtomMapper(atom_map_t *map);

float getRadius(atom_map_t *map, char *atom_name, char *residue_name);
*/
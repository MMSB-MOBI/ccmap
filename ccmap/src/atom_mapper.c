#include "atom_mapper.h"

void create_buffers(char ***names_buffer/*, char ***residue_names_buffer*/, float **radii_buffer){
    int i;

    fprintf(stderr, "Created buffer of length %d\n", MAX_ATOM_PAYLOAD);
    *names_buffer         = malloc( MAX_ATOM_PAYLOAD * sizeof(char**) );
    //*residue_names_buffer = malloc( MAX_ATOM_PAYLOAD * sizeof(char**) );
    for(i = 0 ; i < MAX_ATOM_PAYLOAD ; i++) {
        (*names_buffer)[i]         = malloc( 81 * sizeof(char) );
        //(*residue_names_buffer)[i] = malloc( 81 * sizeof(char) );
    }
    *radii_buffer         = malloc ( MAX_ATOM_PAYLOAD * sizeof(float) );
    fprintf(stderr, "Buffer allocated\n");
}
void destroy_buffers(char **names_buffer, /*char **residue_names_buffer,*/ float *radii_buffer){
    int i;
    fprintf(stderr, "destroying buffer...\n");
    for (i = 0 ; i < MAX_ATOM_PAYLOAD ; i++) {
        free(names_buffer[i]);
        //free(residue_names_buffer[i]);
    }
    free(names_buffer);
    //free(residue_names_buffer);
    free(radii_buffer);
    fprintf(stderr, "buffer free'd\n");
}

atom_map_t *readAtomMapperFromFile(char *filePath, float probeRadius) {
    fprintf(stderr, "readAtomMapperFromFile\n");
    fprintf(stderr, "==>%s\n", filePath);
    FILE *fp;
    atom_map_t *atom_map = NULL;
    int   cur_count = 0;
    float cur_radius;
    char  cur_name[81];
    char  cur_resname[81];
    char prev_resname[81];
    char **atom_name_buffer = NULL;
   // char **res_name_buffer  = NULL;
    float *radii_buffer     = NULL;
   
    cur_name[0] = '\0';
    if ( !(fp = fopen(filePath, "r")) ) {
        fprintf(stderr, "File open error at %s\n", filePath);
        return NULL;
    }

    create_buffers(&atom_name_buffer/*, &res_name_buffer*/, &radii_buffer); 
    atom_map = createAtomMapper();
    while (fscanf(fp, "%s\t%s\t%f\n", cur_resname, cur_name, &cur_radius) != EOF) {
        if( cur_count > 0 &&\
            ! (strcmp(cur_resname, prev_resname) == 0 ) ) {
                addMapGroup(atom_map, atom_name_buffer, prev_resname, radii_buffer, cur_count);    
           // fprintf(stderr, "Changing residue definition\n");
            cur_count = 0;
        }
        strcpy(prev_resname, cur_resname);
        strcpy(atom_name_buffer[cur_count], cur_name);
        radii_buffer[cur_count] = cur_radius;
        cur_count++;
    }
    addMapGroup(atom_map, atom_name_buffer, prev_resname, radii_buffer, cur_count); 

    fclose(fp);
    destroy_buffers(atom_name_buffer,/* res_name_buffer,*/ radii_buffer);
    return atom_map;
}

atom_map_t *destroyAtomMapper(atom_map_t *aMap){
    int i, j;

    for (i = 0 ; i < MAX_RESIDUE_TYPES ; i++) {
        free(aMap->map[i].atom_payload_list);
        for(j = 0 ; j < aMap->map[i].names_length ; j++)
            free(aMap->map[i].names[j]);
        free(aMap->map[i].names);
    }
    free(aMap);

    return NULL;
}

atom_map_t *createAtomMapper() {

    atom_map_t *aMap = malloc(sizeof(atom_map_t));
    aMap->length      = 0;
    aMap->probeRadius = DEFAULT_RADIUS;
    aMap->map         = malloc(MAX_RESIDUE_TYPES * sizeof(atom_payload_map_t));

    return aMap;
}

void addMapGroup(atom_map_t *aMap, char **atom_names, char *residue_names, float *radii, int nb_elem) {
    fprintf(stderr, "addMapGroup for following %s [%d nb_elem]:\n", residue_names, nb_elem);
    int c = 0;
    for (c = 0 ; c < nb_elem ; c++ )
        fprintf(stderr, "%s\t%f\n", atom_names[c], radii[c]);


    atom_payload_map_t *newGroup = &(aMap->map[aMap->length]);
    newGroup->atom_payload_length = nb_elem;
    newGroup->atom_payload_list = malloc(nb_elem * sizeof(atom_payload_t));

    int i = 0;
    int nb_resname = 1;
    char names_buffer[81];

    while(residue_names[i] != '\0') {
        if(residue_names[i] == '|')
            nb_resname++;
        i++;
    }
    newGroup->names_length = nb_resname;
    newGroup->names = malloc(nb_resname * sizeof(char*));
    for (i = 0 ; i < nb_resname ; i++)
        newGroup->names[i] = malloc(81 * sizeof(char));
    
    i = 0;
    nb_resname = 0;
    int offset = 0;
    int j;
    while(residue_names[i] != '\0') {
        
        if(residue_names[i] == '|') {
            for ( j = offset ; offset + j < i ; j++)
                names_buffer[j - offset] = residue_names[offset + j];
            names_buffer[i - offset] = '\0';
            strcpy(newGroup->names[nb_resname], names_buffer);
            offset = i + 1;        
            nb_resname++;
        }
        i++;
    }

    int iAtom;
    for (iAtom = 0 ; iAtom < nb_elem ; iAtom++) {
        newGroup->atom_payload_list[iAtom].radius = radii[iAtom];
        strcpy(newGroup->atom_payload_list[iAtom].name, atom_names[iAtom]);
        
    }
    aMap->length++;
}

float getRadius(atom_map_t *aMap, char *atom_name, char *residue_name) {
    fprintf("getRadius for %s %s\n", residue_name, atom_name);

    int iRes, iAtom, iAltName;
    for (iRes = 0 ; iRes < aMap->length ; iRes++) {
        for(iAltName = 0 ; iAltName <= aMap->map[iRes].names_length ; iAltName++) {
            if( strcmp( aMap->map[iRes].names[iAltName], residue_name ) == 0) {
                for(iAtom = 0 ; iAtom <  aMap->map[iRes].atom_payload_length ; iAtom++)
                    if ( strcmp(aMap->map[iRes].atom_payload_list[iAtom].name, atom_name) == 0 )
                        return aMap->map[iRes].atom_payload_list[iAtom].radius;

                }
            fprintf(stderr, "Unknown atom  %s within residue %s\n", atom_name, residue_name);
            return 0.0;
        }
    }
    fprintf(stderr, "Unknown residue name  %s %s\n", residue_name, atom_name);
    return 0.0;

}

/*
atom_map_t *readAtomMapperFromFile(char *filePath, float probeRadius);
void readMapperInputs(FILE *fp, char **atom_names, char **residue_names, float **radii, float probeRadius);
atom_map_t *createAtomMapper(char **atom_names, char **residue_names, float **radii);
atom_map_t *destroyAtomMapper(atom_map_t *map);

float getRadius(atom_map_t *map, char *atom_name, char *residue_name);
*/
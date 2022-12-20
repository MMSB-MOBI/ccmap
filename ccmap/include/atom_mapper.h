//File: atom_mapper.h
#ifndef ATOM_MAPPER_H
#define ATOM_MAPPER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "miscellaneous.h"

#define MAX_ATOM_PAYLOAD 50 // Max number ot atom beads/atoms per residue
#define MAX_RESIDUE_TYPES 1024 // Max number of residue/atom groups
#define DEFAULT_RADIUS 1.4 // default probe radius
typedef struct {
    char name[81];
    float radius;
} atom_payload_t;

typedef struct {
    char **names; // Residue names and alt names
    uint16_t names_length;
    uint16_t atom_payload_length;
    atom_payload_t *atom_payload_list;
}
atom_payload_map_t;

typedef struct {
    float maxRadius;
    uint16_t length;
    atom_payload_map_t *map;
    float probeRadius;
} atom_map_t;


atom_map_t *readAtomMapperFromFile(char *filePath, float probeRadius);
atom_map_t *destroyAtomMapper(atom_map_t *map);
void create_buffers(char ***names_buffer/*, char ***residue_names_buffer*/, float **radii_buffer);
void destroy_buffers(char **names_buffer/*, char ***residue_names_buffer*/, float *radii_buffer);

atom_map_t *createAtomMapper();
void addMapGroup(atom_map_t *map, char **atom_names, char *residue_name, float *radii, int nb_elem);
/*
void readMapperInputs(FILE *fp, char **atom_names, char **residue_names, float **radii, float probeRadius);
atom_map_t *createAtomMapper(char **atom_names, char **residue_names, float **radii);
*/

float getRadius(atom_map_t *map, char *atom_name, char *residue_name);
void atomMapperPrint(atom_map_t *aMap);
atom_payload_map_t *getPayloadMap(atom_map_t * aMap, char *resname);

#endif

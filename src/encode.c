#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mesh.h"
#ifndef ENCODE
#include "encode.h"
#define ENCODE
#endif
#ifdef AS_PYTHON_EXTENSION
#include <Python.h>
#endif

unsigned int *copyTable(unsigned int *table, int lenTable ){
  unsigned int *newTable = NULL;
  newTable      = malloc( lenTable * sizeof(int) );
  if (newTable != NULL){
    for (int i = 0; i < lenTable; i++){
      newTable[i] = table[i];
    }
  }
  return newTable;
}

/*
By convention the number of columns is taken from the jResidueList if it is not NULL
*/
unsigned int *encodeContactMap(residue_t *iResidueList,         \
                               residue_t *jResidueList, unsigned int *totalContacts
){
    // Initiate table with maximal size
unsigned int *table    = NULL;
unsigned int *newTable = NULL;

residue_t *srcResidue = NULL;
residue_t *tgtResidue = NULL;
size_t tableCapacity = TABLE_CHUNCK_SZ;

unsigned int jLen = jResidueList != NULL ?\
                    residueListLen(jResidueList) :\
                    residueListLen(iResidueList) ;

table = malloc( tableCapacity * sizeof(int));
*totalContacts = 0;
#ifdef DEBUG
  FILE *fp = fopen("encodeContactMap.log", "w"); 

char res1[81];
char res2[81];
#endif

srcResidue = iResidueList;
while (srcResidue != NULL){
  #ifdef DEBUG
    stringifyResidue(srcResidue, res1);
    fprintf(fp, "Looking at %s\n", res1);
  #endif
    if (srcResidue->nContacts > 0) {
        for (int i = 0; i < srcResidue->nContacts; i++){
          if (*totalContacts == (unsigned int)tableCapacity) {
            tableCapacity += TABLE_CHUNCK_SZ;
            table = realloc(table, (size_t)tableCapacity * sizeof (unsigned int) );
          }
          tgtResidue = srcResidue->contactResidueList[i];
          table[*totalContacts] = ENCODE_IJ2K(srcResidue->index, tgtResidue->index, jLen);
          *totalContacts += 1;
          #ifdef DEBUG
          stringifyResidue(tgtResidue, res2);
          fprintf(fp, "%s -- %s\n", res1, res2);
          fprintf(fp, "## %d %d (%d) => %d\n", srcResidue->index, tgtResidue->index, jLen, ENCODE_IJ2K(srcResidue->index, tgtResidue->index, jLen));
          #endif
        }
    }
    srcResidue = srcResidue->nextResidueList;
}

if(jResidueList != NULL) {
  // Now we browse through the ligand, we swap the call to ENCODE
    srcResidue = jResidueList;
    while (srcResidue != NULL){
        #ifdef DEBUG
        stringifyResidue(srcResidue, res2);
        fprintf(fp, "Looking at %s\n", res2);
        #endif
        if (srcResidue->nContacts > 0) {
            for (int i = 0; i < srcResidue->nContacts; i++){
              if (*totalContacts == (unsigned int)tableCapacity) {
                tableCapacity += TABLE_CHUNCK_SZ;
                table = realloc(table, (size_t)tableCapacity * sizeof (unsigned int) );
              }

              tgtResidue = srcResidue->contactResidueList[i];
              table[*totalContacts] = (unsigned int)ENCODE_IJ2K(tgtResidue->index, srcResidue->index, jLen);
              *totalContacts += 1;
             
              #ifdef DEBUG
              stringifyResidue(tgtResidue, res1);
              fprintf(fp, "%s -- %s\n", res1, res2);
              fprintf(fp, "VS## %d %d (%d) => %d\n", tgtResidue->index, srcResidue->index, jLen,\
                                        ENCODE_IJ2K(tgtResidue->index, srcResidue->index, jLen));
              #endif
            }
        }
        srcResidue = srcResidue->nextResidueList;
    }
  }
#ifdef DEBUG
  fclose(fp);
#endif
// Resize table with copyTable and free original table
newTable = copyTable(table, *totalContacts);
free(table);
#ifdef DEBUG
fprintf(stderr, "TT CC %d\n", *totalContacts);
#endif
return newTable;
}

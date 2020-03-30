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

void printTable(int *ContactList, unsigned int len){
  printf("Contact Table :");
  for (unsigned int i=0; i<len; i++){
  printf(" %d ;" ,ContactList[i]);
  }
  printf("\n");
}


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

unsigned int iLen = residueListLen(iResidueList);
unsigned int jLen = jResidueList != NULL ?\
                    residueListLen(jResidueList) :\
                    residueListLen(iResidueList) ;

table = malloc( (iLen * jLen) * sizeof(int));
*totalContacts = 0;

srcResidue = iResidueList;
while (srcResidue->nextResidueList!= NULL){
    if (srcResidue->nContacts > 0) {
        for (int i = 0; i < srcResidue->nContacts; i++){
          tgtResidue = srcResidue->contactResidueList[i];
          table[*totalContacts] = ENCODE_IJ2K(srcResidue->index, tgtResidue->index, jLen);
          *totalContacts += 1;
        }
    }
    srcResidue = srcResidue->nextResidueList;
}

if(jResidueList != NULL) {
  // Now we browse through the ligand, we swap the call to ENCODE
    srcResidue = jResidueList;
    while (srcResidue->nextResidueList != NULL){
        if (srcResidue->nContacts > 0) {
            for (int i = 0; i < srcResidue->nContacts; i++){
              tgtResidue = srcResidue->contactResidueList[i];
              table[*totalContacts] = (unsigned int)ENCODE_IJ2K(tgtResidue->index, srcResidue->index, jLen);
              *totalContacts += 1;
            }
        }
        srcResidue = srcResidue->nextResidueList;
    }
  }

// Resize table with copyTable and free original table
newTable = copyTable(table, *totalContacts);
free(table);
return newTable;
}
/*
int *encodeContactMap(residue_t *ResidueList, int lenLigList, int lenRecList, unsigned int *finalLen){
    // Initiate table with maximal size
    int *table    = NULL;
    int *newTable = NULL;
    residue_t *next_residue=NULL;
    residue_t *residue= ResidueList;
    residue_t *contact= NULL;
    
    table = malloc( (lenLigList*lenRecList) * sizeof(int));
    int o = -1;
    *finalLen = 0;
    while (residue->nextResidueList!= NULL){
      // Check Residues indexes increment, else there is a change of molecule
        next_residue=residue->nextResidueList;
        if (next_residue->index < residue->index){break;}

      // Check there are contacts to store
        int nContacts= residue->nContacts;
        if(residue->nContacts > 0) {
            for (int i=0; i<nContacts; i++){
              o++;
              contact = residue->contactResidueList[i];
             // table[o]=contactIndex(residue->index,contact->index, lenLigList);
              table[o] = ENCODE_IJ2K(residue->index,contact->index, lenLigList);
              // printf(" %d . Contact : ind1 = %d, ind2 = %d  , => %d \n", o, residue->index, contact->index, table[o]);
            }
        }
        residue=next_residue;
    }
    *finalLen=o+1;
    // Resize table with copyTable and free original table
    newTable= copyTable(table,*finalLen);
    free(table);
    return newTable;
}
*/

int contactIndex(int index1, int index2, int max2){
  return index1 *max2 + index2 ;
}

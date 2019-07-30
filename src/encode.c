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


int chainLen(residue_t *ResidueList){
  int nResidues=1;
  residue_t *residue= ResidueList;
  while (residue->nextResidueList!= NULL){
    nResidues++;
    residue=residue->nextResidueList;
  }
  return nResidues;
}

void printTable(int *ContactList, unsigned int len){
  printf("Contact Table :");
  for (unsigned int i=0; i<len; i++){
  printf(" %d ;" ,ContactList[i]);
  }
  printf("\n");
}


int *copyTable(int *table, int lenTable ){
  int *newTable=NULL;
  newTable= malloc(lenTable*sizeof(int));
  if (newTable!=NULL){
    for (int i=0; i<lenTable; i++){
      newTable[i]=table[i];
    }
  }
  return newTable;
}

int *encodeContactMap(residue_t *ResidueList, int lenLigList, int lenRecList, unsigned int *finalLen){
    // Initiate table with maximal size
    int *table=NULL;
    table=malloc(lenLigList*lenRecList*sizeof(int));
    int *newTable=NULL;
    int o=-1;
    residue_t *next_residue=NULL;
    residue_t *residue= ResidueList;
    residue_t *contact= NULL;
    *finalLen=0;
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
              table[o]=contactIndex(residue->index,contact->index, lenLigList);
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

int contactIndex(int index1, int index2, int max2){
  return index1 *max2 + index2 ;
}
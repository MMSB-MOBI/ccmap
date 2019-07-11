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

void printChain(chainedInt_t *ContactList){
  chainedInt_t *contact= ContactList;
  printf("Contact %d ; " ,contact->index);
  while (contact->nextInt!= NULL){
    contact=contact->nextInt;
    printf("%d ; " ,contact->index);
  }

}
chainedInt_t *encodeContactMap(residue_t *ResidueList, int lenLigList){
    residue_t *next_residue=NULL;
    residue_t *residue= ResidueList;
    residue_t *contact= NULL;
    chainedInt_t *contactIndexList=NULL;
    chainedInt_t *previous_Index = NULL;
    while (residue->nextResidueList!= NULL){


      // Check Residues indexes increment, else there is a change of protein
        next_residue=residue->nextResidueList;
        if (next_residue->index < residue->index){break;}

      // Check there are contacts to store
        int nContacts= residue->nContacts;
        if(residue->nContacts > 0) {
            for (int i=0; i<nContacts; i++){

              contact = residue->contactResidueList[i];
              printf(" Contact : ind1 = %d, ind2 = %d  , => %d \n",residue->index, contact->index, contactIndex(residue->index,contact->index, lenLigList));
              chainedInt_t *Index = NULL;

              // Create the ChainedInt object and allocate memory for it
              Index=createChained_Int(contactIndex(residue->index,contact->index, lenLigList));

              // If there is a previous Index object, set it's nextInt value to current pointer
              printf("Adress %p", Index);
              if (previous_Index!=NULL){previous_Index->nextInt=Index;}

              // Set previous_Index value to current pointer
              previous_Index=Index;

              // Store first pointer
              if (contactIndexList==NULL){ contactIndexList=Index; }

              printf("My value is %d \n ",Index->index);

            }
        }
        residue=next_residue;
    }
    printChain(contactIndexList);
    return(contactIndexList);
}

chainedInt_t *createChained_Int(int n) {
    chainedInt_t *indexInt = NULL;
    indexInt=malloc(sizeof(chainedInt_t));
    if (indexInt!=NULL){
      indexInt->index = n;
      indexInt->nextInt = NULL;
    }
    else{
      exit(2);
    }
    return indexInt;
}
int contactIndex(int index1, int index2, int max2){
  return index1 *max2 + index2 ;
}

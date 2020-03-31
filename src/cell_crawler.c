#include "cell_crawler.h"

/*
Create a cellCrawler struct to perform pariwse cell operations featuring:
- Dual OR not-dual atom pair enumeration 
- Residue OR atom level contact registering 
*/
cellCrawler_t *createCellCrawler(bool atomic, bool dual, double dist) {
    #ifdef DEBUG
    fprintf(stderr, "Starting createCellCrawler(%s, %s, %g)\n", \
                                            atomic ? "true" : "false",\
                                            dual   ? "true" : "false", dist);
    #endif
    cellCrawler_t *cellCrawler = malloc(sizeof(cellCrawler_t));
    cellCrawler->atomPairProcess = &processPairwiseDistance;
    cellCrawler->threshold = dist;
    
    cellCrawler->enumerator = &pairwiseCellEnumerate;
    if (dual)
        cellCrawler->enumerator = &pairwiseCellEnumerateDual;
    cellCrawler->dual = dual;

    updaterStruct_t *updater = malloc(sizeof(updaterStruct_t));
    updater->maxSize = ATOM_CC_LIST_CHUNCK;
    updater->atomContactList = malloc( (size_t)ATOM_CC_LIST_CHUNCK * sizeof(atomPair_t) );
    updater->totalByAtom = 0;
    updater->totalByResidue = 0;
    updater->updaterFn =  &updateResidueContact;
    if (atomic)
        updater->updaterFn  = &updateAtomContact;
    
    cellCrawler->updater = updater;
    return cellCrawler;
}

cellCrawler_t *destroyCellCrawler(cellCrawler_t *cellCrawler) {
    free (cellCrawler->updater->atomContactList);
    free (cellCrawler->updater);
    free (cellCrawler);
    return cellCrawler;
}

void extendCellCrawler(cellCrawler_t *cellCrawler) {
    updaterStruct_t *updater = cellCrawler->updater;
    updater->maxSize += ATOM_CC_LIST_CHUNCK;
    atomPair_t *atomContactListRealloc = realloc( updater->atomContactList, updater->maxSize * sizeof(atomPair_t) );
    assert(atomContactListRealloc != NULL);
    updater->atomContactList = atomContactListRealloc;
}

bool processPairwiseDistance(cellCrawler_t* cellCrawler, atom_t* iAtom, atom_t* jAtom) {
    #ifdef DEBUG
    char iAtomString[81];
    char jAtomString[81];
    stringifyAtom(iAtom, iAtomString);
    stringifyAtom(jAtom, jAtomString);
    fprintf(stderr, "Starting processPairwiseDistance T:%g\n", cellCrawler->threshold);
    #endif
    double currDist = distance(iAtom, jAtom);
    #ifdef DEBUG    
    fprintf(stderr, "DD: %g %g || %s || %s \n", currDist, cellCrawler->threshold, iAtomString, jAtomString);
    #endif
    
    // HERE 
    
    char iAtomString[81];
    char jAtomString[81];
    stringifyAtom(iAtom, iAtomString);
    stringifyAtom(jAtom, jAtomString);
    FILE *fp = fopen("contacts_atomic_add.lst", "a");
    fprintf(fp, "DD: %g %g || %s || %s \n", currDist, cellCrawler->threshold, iAtomString, jAtomString);
    

    bool isNewContact = false;
    if(currDist <= cellCrawler->threshold) {     
        #ifdef DEBUG
        fprintf(stderr, "Contact found at %g A\n", currDist);
        #endif
        isNewContact = cellCrawler->updater->updaterFn(cellCrawler, iAtom, jAtom, currDist);
        if(isNewContact) {// Usefull to resize atomContactList only, ...
            // HERE
            fprintf(fp, "%s CC_DIST %s  %s  ==> %.2g\n", cellCrawler->dual ? "dual" : "not_dual", iAtomString, jAtomString, currDist);
            /*
            FILE *fp = fopen("contacts.lst", "a");
            fprintf(fp, "%s CC_DIST %s  %s  ==> %.2g\n", cellCrawler->dual ? "dual" : "not_dual", iAtomString, jAtomString, currDist);
            fclose(fp);
        */
            // THERE
            //cellCrawler->updater->totalByAtom++;
            #ifdef DEBUG 
                fprintf(stderr, "%s CC_DIST %s  %s  ==> %.2g\n", cellCrawler->dual ? "dual" : "", iAtomString, jAtomString, currDist);
                #ifdef AS_PYTHON_EXTENSION
                    PySys_WriteStdout("%s [[Dnum %d]] %s  %s  ==> %.2g\n", cellCrawler->dual ? "dual" : "", iAtomString, jAtomString, currDist);
                #endif
            #endif
        }
    }
    #ifdef DEBUG 
    fprintf(stderr, "Exiting processPairwiseDistance\n");
    #endif

    fclose(fp);

    return isNewContact;
}


bool updateAtomContact(cellCrawler_t *cellCrawler, atom_t *a, atom_t *b, double dist) {
    #ifdef DEBUG 
    fprintf(stderr, "Entering updateAtomContact\n");
    #endif
    updaterStruct_t *updater = cellCrawler->updater;
    if(updater->totalByAtom == updater->maxSize)
        extendCellCrawler(cellCrawler);

    updater->atomContactList[updater->totalByAtom].a    = a;
    updater->atomContactList[updater->totalByAtom].b    = b;
    updater->atomContactList[updater->totalByAtom].dist = dist;
    #ifdef DEBUG 
    fprintf(stderr, "Exiting updateAtomContact\n");
    #endif
    updater->totalByAtom++; // For mem managment
    return true;
}

bool updateResidueContact(cellCrawler_t *cellCrawler, atom_t *iAtom, atom_t *jAtom, double dist) {
    #ifdef DEBUG 
    fprintf(stderr, "Starting updateResidueContact\n");
    #endif
    if(iAtom->belongsTo == jAtom->belongsTo) return false;

    residue_t *iResidue = iAtom->belongsTo;
    residue_t *jResidue = jAtom->belongsTo;

   /*
    FILE *fp = fopen("contacts_residue_add.lst", "a");
    char res1[81];
    char res2[81];
    stringifyResidue(iResidue, res1);
    stringifyResidue(jResidue, res2);
    */
    FILE *fp = fopen("contacts_residue_add.lst", "a");
    char iAtomString[81];
    char jAtomString[81];
    stringifyAtom(iAtom, iAtomString);
    stringifyAtom(jAtom, jAtomString);
    char res1[81];
    stringifyResidue(iResidue, res1);
    char res2[81];
    stringifyResidue(jResidue, res2);
    fprintf(stderr, "%s :: %s -URC- %s :: %s\n", iAtomString, res1, res2, jAtomString);
    fprintf(fp, "%s :: %s -URC- %s :: %s\n", iAtomString, res1, res2, jAtomString);
   
    if (!cellCrawler->dual) { // We order only if its within same pdbrecord
        iResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? iAtom->belongsTo : jAtom->belongsTo;
        jResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? jAtom->belongsTo : iAtom->belongsTo;
    }
    for (int i = 0; i < iResidue->nContacts; i++) {
        if(iResidue->contactResidueList[i] == jResidue) {
            fclose(fp);
          //  fprintf(fp, "CONTACT ALREADY KNOWN between residues %s -- %s\n", res1, res2);
           // fprintf(fp, "Contact already knwow between residues indexed %d,%d\n", iResidue->index, jResidue->index);
            return false;
        }
    }
    for (int j = 0; j < jResidue->nContacts; j++) {
        if(jResidue->contactResidueList[j] == iResidue) {
            fclose(fp);
          //  fprintf(fp, "CONTACT ALREADY KNOWN between residues %s -- %s\n", res1, res2);
           // fprintf(fp, "Contact already knwow between residues indexed %d,%d\n", iResidue->index, jResidue->index);
            return false;
        }
    }
    //HERE
    
    fprintf(fp, "ADDING a new contact between residues %s -- %s\n", res1, res2);
    fclose(fp);
    //THERE
    
    iResidue->nContacts++;
    iResidue->contactResidueList = realloc( iResidue->contactResidueList, iResidue->nContacts * sizeof(residue_t*) );
    iResidue->contactResidueList[iResidue->nContacts - 1] = jResidue;
    #ifdef DEBUG
    fprintf(stderr, "Exiting updateResidueContact\n");
    #endif
    cellCrawler->updater->totalByResidue++; // For logging
    return true;
}

/*
    PairwiseCellEnumerator functions have following prototype
    void (*cellEnumerator)(cellCrawler_t *, 
                           cell_t *, cell_t *);
    pairwiseCellEnumerate
    pairwiseCellEnumerate_DUAL
*/
void pairwiseCellEnumerate(cellCrawler_t *cellCrawler, cell_t *refCell, cell_t *targetCell) {    
    #ifdef DEBUG                        
    fprintf(stderr, "Starting pairwiseCellEnumerate\n");
    fprintf(stderr, "\n*********\nPairwise cell atom enumeration: [%d %d %d] %d elts // [%d %d %d] %d elts\n", \
                                                        refCell->i, refCell->j, refCell->k, refCell->memberCount,\
                                                        targetCell->i, targetCell->j, targetCell->k, targetCell->memberCount);
    #ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("\n*********\nPairwise cell atom enumeration: [%d %d %d] %d elts // [%d %d %d] %delts\n", \
                                                        refCell->i, refCell->j, refCell->k, refCell->memberCount,\
                                                        targetCell->i, targetCell->j, targetCell->k, targetCell->memberCount);
   #endif
    char iAtomString[81];
    char jAtomString[81];
    #endif
    
    if(refCell->memberCount == 0 || targetCell->memberCount == 0) return;

    atom_t *iAtom, *jAtom;
    iAtom = refCell->members;

    while(iAtom != NULL) {

    #ifdef DEBUG
        stringifyAtom(iAtom, iAtomString);
    #endif
        jAtom = targetCell->members;
        while(jAtom != NULL) {
            if (jAtom != iAtom) {
                #ifdef DEBUG
                stringifyAtom(jAtom, jAtomString);
                fprintf(stderr, "calling atomPairProcess for  %s VS %s\n", iAtomString, jAtomString);
                #endif
                cellCrawler->atomPairProcess(cellCrawler, iAtom, jAtom);
            }
            jAtom = jAtom->nextCellAtom;
        }
        iAtom = iAtom->nextCellAtom;
    }
}

void pairwiseCellEnumerateDual(cellCrawler_t *cellCrawler, cell_t *refCell, cell_t *targetCell) {
    atom_t *iAtom, *jAtom;
#ifdef DEBUG
    char iAtomString[81];
    char jAtomString[81];
#endif
    if(refCell->memberCount == 0 || targetCell->memberCount == 0) return;
#ifdef DEBUG
    printf("\n*********\nDUAL Pairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("DUAL Pairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#endif
#endif


    iAtom = refCell->iMembers;
    while(iAtom != NULL) {

    #ifdef DEBUG
        stringifyAtom(iAtom, iAtomString);
    #endif
        jAtom = targetCell->jMembers;
        while(jAtom != NULL) {
            if (jAtom != iAtom) {
                #ifdef DEBUG
                    stringifyAtom(jAtom, jAtomString);
                #endif
                cellCrawler->atomPairProcess(cellCrawler, iAtom, jAtom);                
            }
            jAtom = jAtom->nextCellAtom;
        }
        iAtom = iAtom->nextCellAtom;
    }
// Reverse i/j members lookup
    iAtom = refCell->jMembers;
    while(iAtom != NULL) {

    #ifdef DEBUG
        stringifyAtom(iAtom, iAtomString);
    #endif
        jAtom = targetCell->iMembers;
        while(jAtom != NULL) {
            if (jAtom != iAtom) {
                #ifdef DEBUG
                    stringifyAtom(jAtom, jAtomString);
                #endif
                cellCrawler->atomPairProcess(cellCrawler, iAtom, jAtom);
            }
            jAtom = jAtom->nextCellAtom;
        }
        iAtom = iAtom->nextCellAtom;
    }
#ifdef DEBUG
    fprintf(stderr, "Exiting pairwiseCellEnumerate\n");
#endif
}

double distance(atom_t *iAtom, atom_t *jAtom) {
    FILE *fp = fopen("shadow.lst", "a");
    char a1[100];
    char a2[100];
    double _ = sqrt( (iAtom->x - jAtom->x) * (iAtom->x - jAtom->x) + (iAtom->y - jAtom->y) * (iAtom->y - jAtom->y) + (iAtom->z - jAtom->z) * (iAtom->z - jAtom->z));
    stringifyAtom(iAtom, a1);
    stringifyAtom(jAtom, a2);

    fprintf(fp, "%s %s %.5g\n", a1, a2, _);
    fclose(fp);


    return sqrt( (iAtom->x - jAtom->x) * (iAtom->x - jAtom->x) + (iAtom->y - jAtom->y) * (iAtom->y - jAtom->y) + (iAtom->z - jAtom->z) * (iAtom->z - jAtom->z) );
}

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

// iAtom is guaranted to be from 1st body coordinate sets, or single
// jAtom is guaranted to be from 2 body coordinate sets, or single
bool processPairwiseDistance(cellCrawler_t* cellCrawler, atom_t* iAtom, atom_t* jAtom) {
    double currDist = distance(iAtom, jAtom);
    
    #ifdef DEBUG    
    fprintf(stderr, "Starting processPairwiseDistance T:%g\n", cellCrawler->threshold);
    char iAtomString[81];
    char jAtomString[81];
    stringifyAtom(iAtom, iAtomString);
    stringifyAtom(jAtom, jAtomString);
    FILE *fp = fopen("contacts_atomic_add.lst", "a");
    fprintf(fp, "DD: %g %g || %s || %s \n", currDist, cellCrawler->threshold, iAtomString, jAtomString);
    #endif

    bool isNewContact = false;
    if(currDist <= cellCrawler->threshold) {     
        #ifdef DEBUG
        fprintf(stderr, "Contact found at %g A\n", currDist);
        #endif
        isNewContact = cellCrawler->updater->updaterFn(cellCrawler, iAtom, jAtom, currDist);
        if(isNewContact) {
            #ifdef DEBUG 
                fprintf(fp, "%s CC_DIST %s  %s  ==> %.2g\n", cellCrawler->dual ? "dual" : "not_dual", iAtomString, jAtomString, currDist);
            #endif
        }
    }
    #ifdef DEBUG 
    fprintf(stderr, "Exiting processPairwiseDistance\n");
    fclose(fp);
    #endif

    return isNewContact;
}

// a atom is guaranted to be from 1st body coordinate sets, or single
// b atom is guaranted to be from 2 body coordinate sets, or single
bool updateAtomContact(cellCrawler_t *cellCrawler, atom_t *a, atom_t *b, double dist) {
    #ifdef DEBUG 
    fprintf(stderr, "Entering updateAtomContact\n");
    #endif
    updaterStruct_t *updater = cellCrawler->updater;
    if(updater->totalByAtom == updater->maxSize)
        extendCellCrawler(cellCrawler);
    updater->atomContactList[updater->totalByAtom].next = NULL;
    updater->atomContactList[updater->totalByAtom].a    = a;
    updater->atomContactList[updater->totalByAtom].b    = b;
    updater->atomContactList[updater->totalByAtom].dist = dist;
    if(updater->totalByAtom > 0)
        updater->atomContactList[updater->totalByAtom - 1 ].next = \
        &updater->atomContactList[updater->totalByAtom];
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

  
#ifdef DEBUG 
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
#endif
    if (!cellCrawler->dual) { // We order only if its within same pdbrecord
        iResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? iAtom->belongsTo : jAtom->belongsTo;
        jResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? jAtom->belongsTo : iAtom->belongsTo;
    }
    for (int i = 0; i < iResidue->nContacts; i++) {
        if(iResidue->contactResidueList[i] == jResidue) {
            #ifdef DEBUG 
            fclose(fp);
            #endif
            return false;
        }
    }
    for (int j = 0; j < jResidue->nContacts; j++) {
        if(jResidue->contactResidueList[j] == iResidue) {
            #ifdef DEBUG 
            fclose(fp);
            #endif
            return false;
        }
    }
    #ifdef DEBUG   
    fprintf(fp, "ADDING a new contact between residues %s -- %s\n", res1, res2);
    fclose(fp);
    #endif
    
    iResidue->nContacts++;
    iResidue->contactResidueList = realloc( iResidue->contactResidueList, iResidue->nContacts * sizeof(residue_t*) );
    iResidue->contactResidueList[iResidue->nContacts - 1] = jResidue;
    #ifdef DEBUG
    fprintf(stderr, "Exiting updateResidueContact\n");
    #endif
    cellCrawler->updater->totalByResidue++; // For logging
    return true;
}

void pairwiseCellEnumerate(cellCrawler_t *cellCrawler, cell_t *refCell, cell_t *targetCell) {    
    #ifdef DEBUG                        
    fprintf(stderr, "Starting pairwiseCellEnumerate\n");
    fprintf(stderr, "\n*********\nPairwise cell atom enumeration: [%d %d %d] %d elts // [%d %d %d] %d elts\n", \
                                                        refCell->i, refCell->j, refCell->k, refCell->memberCount,\
                                                        targetCell->i, targetCell->j, targetCell->k, targetCell->memberCount);
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
                cellCrawler->atomPairProcess(cellCrawler, jAtom, iAtom);
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
    #ifdef DEBUG 
    FILE *fp = fopen("shadow.lst", "a");
    char a1[100];
    char a2[100];
    double _ = sqrt( (iAtom->x - jAtom->x) * (iAtom->x - jAtom->x) + (iAtom->y - jAtom->y) * (iAtom->y - jAtom->y) + (iAtom->z - jAtom->z) * (iAtom->z - jAtom->z));
    stringifyAtom(iAtom, a1);
    stringifyAtom(jAtom, a2);

    fprintf(fp, "%s %s %.5g\n", a1, a2, _);
    fclose(fp);
    #endif

    return sqrt( (iAtom->x - jAtom->x) * (iAtom->x - jAtom->x) + (iAtom->y - jAtom->y) * (iAtom->y - jAtom->y) + (iAtom->z - jAtom->z) * (iAtom->z - jAtom->z) );
}

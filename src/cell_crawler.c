#include "cell_crawler.h"

/*
Create a cellCrawler struct to perform pariwse cell operations featuring:
- Dual OR not-dual atom pair enumeration 
- Residue OR atom level contact registering 
*/
cellCrawler_t createCellCrawler(bool atomic, bool dual, double dist) {
    cellCrawler_t cellCrawler;
    cellCrawler.atomPairProcess = &processPairwiseDistance;
    cellCrawler.threshold = dist;
    
    cellCrawler.enumerator = &pairwiseCellEnumerate;
    if (dual)
        cellCrawler.enumerator = &pairwiseCellEnumerateDual;
    cellCrawler.dual = dual;

    updaterStruct_t updater;
    updater.maxSize = ATOM_CC_LIST_CHUNCK;
    updater.atomContactList = malloc( (size_t)ATOM_CC_LIST_CHUNCK * sizeof(atomPair_t) );
    updater.totalByAtom = 0;
    updater.totalByResidue = 0;
    updater.updaterFn =  &updateResidueContact;
    if (atomic)
        updater.updaterFn  = &updateAtomContact;
    
    cellCrawler.updater = &updater;
    return cellCrawler;
}

cellCrawler_t *destroyCellCrawler(cellCrawler_t *cellCrawler) {
    free (cellCrawler->updater->atomContactList);
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
    double currDist = distance(iAtom, jAtom);
    bool isNewContact = false;
    if(currDist <= cellCrawler->threshold) {     
        isNewContact = cellCrawler->updater->updaterFn(cellCrawler, iAtom, jAtom, currDist);
        if(isNewContact) {// Usefull to resize atomContactList only, ...
            cellCrawler->updater->totalByAtom++;
            #ifdef DEBUG
                char iAtomString[81];
                char jAtomString[81];
                stringifyAtom(iAtom, iAtomString);
                stringifyAtom(jAtom, jAtomString);
                printf("%s [[Dnum %d]] %s  %s  ==> %.2g\n", cellCrawler->dual ? "dual" : "", iAtomString, jAtomString, currDist);
                #ifdef AS_PYTHON_EXTENSION
                    PySys_WriteStdout("%s [[Dnum %d]] %s  %s  ==> %.2g\n", cellCrawler->dual ? "dual" : "", iAtomString, jAtomString, currDist);
                #endif
            #endif
        }
    }
    return isNewContact;
}

bool updateAtomContact(cellCrawler_t *cellCrawler, atom_t *a, atom_t *b, double dist) {
    updaterStruct_t *updater = cellCrawler->updater;

    if(updater->totalByAtom == updater->maxSize)
        extendCellCrawler(cellCrawler);

    updater->atomContactList[updater->totalByAtom].a    = a;
    updater->atomContactList[updater->totalByAtom].b    = b;
    updater->atomContactList[updater->totalByAtom].dist = dist;
    
    return true;
}

bool updateResidueContact(cellCrawler_t *cellCrawler, atom_t *iAtom, atom_t *jAtom, double dist) {
    if(iAtom->belongsTo == jAtom->belongsTo) return false;

    residue_t *iResidue = iAtom->belongsTo;
    residue_t *jResidue = jAtom->belongsTo;

    if (!cellCrawler->dual) { // We order only if its within same pdbrecord
        iResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? iAtom->belongsTo : jAtom->belongsTo;
        jResidue = iAtom->belongsTo->index < jAtom->belongsTo->index ? jAtom->belongsTo : iAtom->belongsTo;
    }
    for (int i = 0; i < iResidue->nContacts; i++) {
        if(iResidue->contactResidueList[i] == jResidue) {
          //  printf("Contact already knwow between residues indexed %d,%d\n", iResidue->index, jResidue->index);
            return false;
        }
    }
    //printf("ADDING a new contact between residues indexed %d,%d\n", iResidue->index, jResidue->index);
    iResidue->nContacts++;
    iResidue->contactResidueList = realloc( iResidue->contactResidueList, iResidue->nContacts * sizeof(residue_t*) );
    iResidue->contactResidueList[iResidue->nContacts - 1] = jResidue;
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
    atom_t *iAtom, *jAtom;

    char iAtomString[81];
    char jAtomString[81];
   

    if(refCell->memberCount == 0 || targetCell->memberCount == 0) return;
#ifdef DEBUG
    printf("\n*********\nPairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#ifdef AS_PYTHON_EXTENSION
    PySys_WriteStdout("Pairwise cell atom enumeration: [%d %d %d]// [%d %d %d]\n", refCell->i, refCell->j, refCell->k, targetCell->i, targetCell->j, targetCell->k);
#endif
#endif
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

    char iAtomString[81];
    char jAtomString[81];

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

}

double distance(atom_t *iAtom, atom_t *jAtom) {
    return sqrt( (iAtom->x - jAtom->x) * (iAtom->x - jAtom->x) + (iAtom->y - jAtom->y) * (iAtom->y - jAtom->y) + (iAtom->z - jAtom->z) * (iAtom->z - jAtom->z) );
}

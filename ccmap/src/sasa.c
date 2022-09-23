#include "sasa.h"
#include "mesh.h"

string_t *jsonifySasaResults(sasaResults_t *sasaResults) {
    #ifdef DEBUG
    fprintf(stderr, "Starting jsonifySasaResults over %d sasaResults\n",\
        sasaResults->length);
    #endif
    string_t *jsonString = createString();
    jsonString->append(jsonString, "{freeSASA : [");
    if (sasaResults->length == 0) {
        fprintf(stderr, "sasaResults empty\n");
        jsonString->append(jsonString, "]}");
        return jsonString;
    }
    char buffer[1024];
    residue_sasa_t *residue_sasa = NULL;
    residue_t *residue_curr      = NULL;
    for (int i = 0 ; i < sasaResults->length ; i++) {
       // fprintf(stderr, "-->%d\n", i);
        jsonString->append(jsonString, "{\"residue\":");
        //fprintf(stderr, "toto\n");
        residue_sasa = &sasaResults->residueSasaList[i];
        //fprintf(stderr, "toto2\n");
        residue_curr = residue_sasa->residue;
        //fprintf(stderr, "toto3\n");
        jsonifyResidue(residue_curr, buffer);
        //fprintf(stderr, "toto3a\n");
        jsonString->append(jsonString, buffer);
        //fprintf(stderr, "toto4\n");
        sprintf(buffer, ", \"SASA\": %f, \"norm\": %f, \"frac\": %f]}", \
                       residue_sasa->nominal - residue_sasa->buried,\
                       residue_sasa->nominal, residue_sasa->frac);
        jsonString->append(jsonString, buffer);
        if(i < sasaResults->length - 1)
            jsonString->append(jsonString, ",\n");
    }
    jsonString->append(jsonString, "]}");
   
    return jsonString;
}
/* Compute freeSASA over a list of residues */
sasaResults_t *computeSasaResults(residueList_t *residueList) {
    #ifdef DEBUG
        char residueBuffer[81];           
        fprintf(stderr, "\ncomputeSasaResults:starting\n");
    #endif
    sasaResults_t *sasaResults   = malloc(sizeof(sasaResults_t));
    sasaResults->residueSasaList = malloc( residueList->length *sizeof(residue_sasa_t) );
    float tSurface, bSurface;
    residue_t *currentResidue = residueList->root;
    uint16_t iResidue = 0; 
    #ifdef DEBUG
        fprintf(stderr, "computeSasaResults: Iteration start\n");
    #endif
    residue_sasa_t *currentResidueSasa = NULL;
    while(currentResidue != NULL) {

        currentResidueSasa = &sasaResults->residueSasaList[iResidue];
        currentResidueSasa->residue   = currentResidue;
        currentResidueSasa->nominal   = 0;
        currentResidueSasa->buried    = 0;
    
        #ifdef DEBUG
            fprintf(stderr, "computeSasaResults: Positioned to residue %d over %d residues total\n", \
                            iResidue, residueList->length);
            printResidue(stderr, currentResidue);
        #endif
    
        for (int iAtom = 0 ; iAtom < currentResidue->nAtoms ; iAtom++) {
            #ifdef DEBUG
                fprintf(stderr, "computeSasaResults:%d %d [max is %d]\n", iResidue, iAtom, currentResidue->nAtoms);
            #endif
            computeFiboSphereASA(currentResidue->elements[iAtom].f_grid, &tSurface, &bSurface);
            currentResidueSasa->nominal += tSurface ;
            currentResidueSasa->buried  += bSurface;
        }
        #ifdef DEBUG
            stringifyResidue(currentResidue, residueBuffer);
            fprintf(stderr, "computeSasaResults: assigning %f %f to following residue\n%s\n\n",\
            sasaResults->residueSasaList[iResidue].buried,\
            sasaResults->residueSasaList[iResidue].nominal,\
            residueBuffer);
        #endif
        currentResidueSasa->frac = \
            currentResidueSasa->buried / currentResidueSasa->nominal;

        currentResidue = currentResidue->nextResidueList;
        iResidue++;
    }
    sasaResults->length = iResidue;
#ifdef DEBUG
    fprintf(stderr, "Exiting createSasaResults\n");
#endif

    return sasaResults;
}

sasaResults_t *destroySasaResults(sasaResults_t *sasaResults){
    free(sasaResults->residueSasaList);
    free(sasaResults);

    return sasaResults;
}

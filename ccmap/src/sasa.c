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
    char residue_buffer[81];
    residue_sasa_t *residue_sasa = NULL;
   
    for (int i = 0 ; i < sasaResults->length ; i++) {
        jsonString->append(jsonString, "{\"residue\":");
        residue_sasa = &sasaResults->residueSasaList[i];
        sprintf(residue_buffer,"%c %s %s", residue_sasa->chainID, residue_sasa->resname, residue_sasa->resID);
        jsonString->append(jsonString, residue_buffer);
        
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
        strcpy(currentResidueSasa->resname, currentResidue->resName);
        strcpy(currentResidueSasa->resID, currentResidue->resID);
        currentResidueSasa->chainID  = currentResidue->chainID;
        currentResidueSasa->res_index  = currentResidue->index;
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

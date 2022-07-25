#include "sasa.h"
#include "mesh.h"

/* Computing sasa w/out contact as DVL */
sasaResults_t *sasaCore(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, float probeRadius) {

    residue_t *iResidueList                     = createResidueList(iAtomList);
    residue_t *jResidueList = jAtomList != NULL ? createResidueList(jAtomList) : NULL;

    double step = probeRadius * 2 + VDW_MAX * 2;
    meshContainer_t *meshContainer = createMeshContainer(iAtomList, iAtom, jAtomList, jAtom, step);
   
    cellSasaCrawler_t *cellSasaCrawler = createCellSasaCrawler();
    meshCrawler(meshContainer, NULL, cellSasaCrawler);
    sasaResults_t *results = createSasaResults_t(cellSasaCrawler, iResidueList, jResidueList);
    meshContainer = destroyMeshContainer(meshContainer);
    destroyCellSasaCrawler(cellSasaCrawler);
   
#ifdef DEBUG
    fprintf(stderr, "Exiting sasaCore\n");
#endif

    return results;
}

sasaResults_t *createSasaResults_t(cellSasaCrawler_t *cellSasaCrawler,residue_t *iResidueList, residue_t *jResidueList){

}
sasaResults_t *destroySasaResults_t(sasaResults_t *createSasaResults_t){

}
cellSasaCrawler_t *createCellSasaCrawler() {

}

cellSasaCrawler_t *destroyCellSasaCrawler(cellSasaCrawler_t *createCellSasaCrawler) {


}
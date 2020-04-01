
void ccmap_compute_list_cleanOnExit(ccmapView_t ***ccmapViewList, \
                                    atom_t ***atomListRecList, int **atomListRecSizes, \
                                    atom_t ***atomListLigList, int **atomListLigSizes, \
                                    int nStructPair, bool dual)
{
*atomListRecList      = PyMem_New(atom_t*     , nStructPairs);
*atomListRecSizes     = PyMem_New(size_t*     , nStructPairs);
if(dual) {
    *atomListLigList  = PyMem_New(atom_t*     , nStructPairs);   
    *atomListLigSizes = PyMem_New(size_t*     , nStructPairs);
}
*ccmapViewList        = PyMem_New(ccmapView_t*, nStructPairs);

}

void ccmap_compute_list_cleanOnExit(ccmapView_t **ccmapViewList, \
                                    atom_t **atomListRecList, int *atomListRecSize,\
                                    atom_t **atomListLigList, int *atomListLigSize,\
                                    int nStructPair, bool dual)
{
for (int i = 0 ; i < nStructPair ; i++) {
    destroyAtomList( atomListRecList[i], (int)atomListRecSize[i] );
    if(dual)
        destroyAtomList( atomListLigList[i], (int)atomListLigSize[i] );
}
    
PyMem_Free(atomListRecList);
PyMem_Free(atomListRecSize);
if(dual){ 
    PyMem_Free(atomListLigList);
    PyMem_Free(atomListLigSize);
}
PyMem_Free(ccmapViewList);
}
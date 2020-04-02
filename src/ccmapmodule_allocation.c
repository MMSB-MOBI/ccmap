#include "Python.h"
#include "ccmapmodule_allocation.h"

void ccmap_compute_list_allocate(ccmapView_t ***ccmapViewList, \
                                    atom_t ***atomListRecList, int **atomListRecSizes, \
                                    atom_t ***atomListLigList, int **atomListLigSizes, \
                                    int nStructPair, bool dual)
{
*atomListRecList      = PyMem_New(atom_t*     , nStructPair);
*atomListRecSizes     = PyMem_New(int     , nStructPair);
if(dual) {
    *atomListLigList  = PyMem_New(atom_t*     , nStructPair);   
    *atomListLigSizes = PyMem_New(int     , nStructPair);
}
*ccmapViewList        = PyMem_New(ccmapView_t*, nStructPair);

}

void ccmap_compute_list_cleanOnExit(ccmapView_t **ccmapViewList, \
                                    atom_t **atomListRecList, int *atomListRecSize,\
                                    atom_t **atomListLigList, int *atomListLigSize,\
                                    int nStructPair, bool dual)
{
for (int i = 0 ; i < nStructPair ; i++) {
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

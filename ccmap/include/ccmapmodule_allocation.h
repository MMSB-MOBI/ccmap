#ifndef CCMAP_ALLOC_H_INCLUDED
#define CCMAP_ALLOC_H_INCLUDED
#include <Python.h>
#include <stdlib.h>
#include "mesh.h"
#include "molecular_object.h"
void ccmap_compute_list_allocate(ccmapView_t ***, atom_t ***,  int **, atom_t ***,  int **, int, bool);
void ccmap_compute_list_cleanOnExit(ccmapView_t **, atom_t **, int *, atom_t **, int *, int, bool);
#endif

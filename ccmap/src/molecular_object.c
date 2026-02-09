#include "molecular_object.h"

#ifdef AS_PYTHON_EXTENSION
    #include <Python.h>
#endif

void printResidueList(FILE *stream, residueList_t *residueList) {
    residue_t *curr_residue = residueList->root;
    if(stream == NULL)
        stream = stdout;
    fprintf(stream, "PRINTING RESIDUE LIST\n");

    printResidue(stream, curr_residue);
    while(curr_residue->nextResidueList != NULL) {
        curr_residue = curr_residue->nextResidueList;
        printResidue(stream, curr_residue);
    }
}
/* Deprecated
unsigned int residueListLen(residue_t *ResidueList){
  int nResidues=1;
  residue_t *residue= ResidueList;
  while (residue->nextResidueList!= NULL){
    nResidues++;
    residue=residue->nextResidueList;
  }
  return nResidues;
}
*/
void printResidue(FILE *stream, residue_t *residue) {
    if(stream == NULL)
        stream = stdout;
    atom_t *curr_atom;
    char atomString[81];
    char residueString[81];
    stringifyResidue(residue, residueString);
    fprintf(stream, "RES: %s (%d)\n", residueString, residue->nAtoms);
    curr_atom = residue->elements;
    while(curr_atom->nextResidueAtom != NULL) {
        stringifyAtom(curr_atom, atomString);
        fprintf(stream, "%s\n", atomString);
        curr_atom = curr_atom->nextAtomList;
    }
    stringifyAtom(curr_atom, atomString);
    fprintf(stream, "%s\n", atomString);
    fprintf(stream, "--------------------\n");
}

void printAtomList(atom_t *atomList, FILE *stream) {
    #ifdef DEBUG
    fprintf(stderr, "Starting printAtomList\n");
    #endif
    atom_t *curr_atom = &atomList[0];
    char atomString[1024];
    while(curr_atom->nextAtomList != NULL) {
        stringifyAtom(curr_atom, atomString);
        fprintf(stream == NULL ? stdout : stream,\
                "%s\n", atomString);
        curr_atom = curr_atom->nextAtomList;
    }
    
    stringifyAtom(curr_atom, atomString);
    printf("%s\n", atomString);
    #ifdef DEBUG
    fprintf(stderr, "Exiting printAtomList\n");
    #endif
}

char *stringifyAtomList(atom_t *atomList) {
    
    atom_t *curr_atom = &atomList[0];
    

    char atomString[1024];
    #ifdef DEBUG
        char dbgBuffer[1024];
        int i = 0;
    #endif
    
    string_t *atomListString = createString();
    while(curr_atom->nextAtomList != NULL) {
        #ifdef DEBUG
            sprintf(dbgBuffer, "stringifyAtomList[%d]: %f %f %f\n", i, curr_atom->x, curr_atom->y, curr_atom->z);
            printOnContextStderr(dbgBuffer);
            i++;
        #endif
        stringifyAtom(curr_atom, atomString);
        atomListString->append(atomListString, atomString);
        atomListString->append(atomListString, "\n");
        curr_atom = curr_atom->nextAtomList;
    }
    #ifdef DEBUG
            sprintf(dbgBuffer, "stringifyAtomList(last): %f %f %f\n", curr_atom->x, curr_atom->y, curr_atom->z);
            printOnContextStderr(dbgBuffer);
    #endif
    stringifyAtom(curr_atom, atomString);
    atomListString->append(atomListString, atomString);
    atomListString->append(atomListString, "\n");

    char *atomListChar       = atomListString->toChar(atomListString);
    destroyString(atomListString);
    return atomListChar;
}

void stringifyAtom(atom_t *atom, char *atomString) {
    sprintf(atomString, "%s %s %s %c %f %f %f %f %f", atom->name, atom->resName, atom->resID, atom->chainID, atom->x, atom->y, atom->z, atom->_VDWradius, atom->_radiusASA);    
}
void stringifyResidue(residue_t *residue, char *residueString) {
    sprintf(residueString, "%s %s %c [i=%d, n=%d]", residue->resName, residue->resID, residue->chainID, residue->index, residue->nAtoms);
}

void jsonifyResidue(residue_t *residue, char *jsonString) {
    sprintf(jsonString, "{\"resID\" : \"%s\", \"chainID\":\"%c\"}", residue->resID, residue->chainID);
}

// fixed length JSON array stringrepresentation
void jsonArrayifyAtom(atom_t *atom, char *atomString, bool bCoordinates) {
    if (bCoordinates)
        sprintf(atomString, "[ \"%s\", \"%s\", \"%s\", \"%c\", %g, %g, %g ]", atom->name, atom->resName, atom->resID, atom->chainID, atom->x, atom->y, atom->z);
    else 
        sprintf(atomString, "[ \"%s\", \"%s\", \"%s\", \"%c\" ]", atom->name, atom->resName, atom->resID, atom->chainID);
}
void jsonifyAtomPair(atomPair_t *atomPair, char *jsonString) {
    bool bCoordinates = false;
    char atomA[1024];
    char atomB[1024];
    jsonArrayifyAtom(atomPair->a, atomA, bCoordinates);
    jsonArrayifyAtom(atomPair->b, atomB, bCoordinates);
    sprintf(jsonString, "[ %s, %s, %5.3g ]", atomA, atomB, atomPair->dist);
}

/*
Replace the x,yz values of the 2nd list by those of the 1st list
Fibo grid managment needs checking, should be ok 
*/
bool applyCoordinates(atom_t *atomListFrom, atom_t *atomListTo) {
    
    while (atomListFrom != NULL && atomListTo != NULL ) {
        float dx, dy, dz;
        dx = atomListFrom->x - atomListTo->x;
        dy = atomListFrom->y - atomListTo->y;
        dz = atomListFrom->z - atomListTo->z;
        
        atomListTo->x = atomListFrom->x;
        atomListTo->y = atomListFrom->y;
        atomListTo->z = atomListFrom->z;

        if(atomListFrom->f_grid != NULL) 
            updateFiboGrid(atomListTo->f_grid, dx, dy, dz);
        
        atomListFrom = atomListFrom->nextAtomList;
        atomListTo   = atomListTo->nextAtomList;
    }
    
    return atomListTo == NULL && atomListFrom == NULL;
}

/*
Update atomList coordinates, translating fiboSphere
*/
void updateCoordinates(atom_t *atomList, coordinates_t *coordinates) {
    
    int curr_index = 0;
    float dx, dy, dz;
    atom_t *currAtom = atomList;
    while (currAtom != NULL) {
        dx = coordinates[curr_index].x - currAtom->x;
        dy = coordinates[curr_index].y - currAtom->y;
        dz = coordinates[curr_index].z - currAtom->z;
        currAtom->x = coordinates[curr_index].x;
        currAtom->y = coordinates[curr_index].y;
        currAtom->z = coordinates[curr_index].z;

        if(currAtom->f_grid != NULL) 
            updateFiboGrid(currAtom->f_grid, dx, dy, dz);
        
        currAtom = currAtom->nextAtomList;
        curr_index++;
    }
}


residue_t *createResidue(atom_t *atom, int n) {
    #ifdef DEBUG
    FILE *fp=fopen("createResidue.log", "a"); 
    char atomString[81];
    stringifyAtom(atom, atomString);
    fprintf(fp, "Starting createResidue with %s\n", atomString);
    fclose(fp);
    #endif
    residue_t *residue = malloc(sizeof(residue_t));
    residue->nAtoms = 1;
    residue->nContacts = 0;
    residue->resName = malloc( (strlen(atom->resName) + 1) * sizeof(char) );
    strcpy(residue->resName, atom->resName);
    residue->resID = malloc( (strlen(atom->resID) + 1) * sizeof(char) );
    strcpy(residue->resID, atom->resID);
    residue->ext_chainID = NULL;
    if (atom->ext_chainID != NULL) {
        residue->ext_chainID = malloc( (strlen(atom->ext_chainID) + 1) * sizeof(char));
        strcpy(residue->ext_chainID, atom->ext_chainID);
    }
    residue->chainID = atom->chainID;
    residue->elements = atom;
    residue->elements->nextResidueAtom = NULL;
    residue->contactResidueList = NULL;
    residue->index = n;
    residue->prevResidueList = NULL;
    residue->nextResidueList = NULL;

    return residue;
}

// operate there GL
residueList_t *createResidueList(atom_t * atomList) {
    #ifdef DEBUG
    fprintf(stderr, "Entering createResidueList\n");
    #endif
    residueList_t *residueList = malloc(sizeof(residueList_t));
    residueList->length = 0; 
    //
       
    residue_t *residue_head = NULL;
    atom_t *curr_atom = atomList;//&atomList[0];
    atom_t *prev_atom = NULL;
    residueList->root   = createResidue(curr_atom,  residueList->length);

    residue_head = residueList->root;
    curr_atom->belongsTo = residue_head;
    prev_atom = curr_atom;
    curr_atom = curr_atom->nextAtomList;

    bool isSameChain;
    while(curr_atom != NULL) {
        isSameChain = true;
        if (curr_atom->chainID != residue_head->chainID)
            isSameChain = false;
        if(curr_atom->ext_chainID != NULL)
            if (strcmp(curr_atom->ext_chainID, residue_head->ext_chainID) !=0)
                isSameChain = false;
        if (strcmp(curr_atom->resID, residue_head->resID) == 0\
            && isSameChain) {
            prev_atom->nextResidueAtom = curr_atom;
            residue_head->nAtoms += 1;
        } else {
            residueList->length++;
            residue_head->nextResidueList = createResidue(curr_atom, residueList->length);
            residue_head->nextResidueList->prevResidueList = residue_head;
            residue_head = residue_head->nextResidueList;
        }
        curr_atom->nextResidueAtom = NULL;
        prev_atom = curr_atom;
        curr_atom->belongsTo = residue_head;

        curr_atom = curr_atom->nextAtomList;
    }
    // GL patches 09/2022 set to actual length, seems to be w/out effect as all 
    // residueList operations seems to employ the chain list pointers
    residueList->length += 1;
    //--
#ifdef DEBUG
    fprintf(stderr, "Created %d residues, exiting createResidueList\n", residueList->length + 1);
#endif
    //
    return residueList;
}

// Fuse the 2nd list into the first, freeing the 2nd one and returning its pointer as NULL
residueList_t *fuseResidueLists(residueList_t *iResidueList, residueList_t *jResidueList) {
    residue_t *head = iResidueList->root;
    while(head->nextResidueList != NULL)
        head = head->nextResidueList;
    head->nextResidueList = jResidueList->root;
    jResidueList->root->prevResidueList = head;
    iResidueList->length += jResidueList->length;
    // we free the wraping jList structure
    free(jResidueList);
    return jResidueList;
}

// Move to end of the chain and free residues in a backward motion
residueList_t *destroyResidueList(residueList_t *residueList) {
#ifdef DEBUG
    fprintf(stderr, "Destroying residue list\n");
#endif
    assert(residueList != NULL);
    residue_t *curr_residue = residueList->root;
    while (curr_residue->nextResidueList != NULL) {
        curr_residue = curr_residue->nextResidueList;
    }
    while (curr_residue->prevResidueList != NULL) {
        curr_residue = curr_residue->prevResidueList;
        curr_residue->nextResidueList = destroyResidue(curr_residue->nextResidueList);
    }
    residueList->root = destroyResidue(curr_residue);
#ifdef DEBUG
    fprintf(stderr, "Destroying residue list: OK\n");
#endif
    free(residueList);
    return NULL;
}

residue_t *destroyResidue(residue_t *residue) {
#ifdef DEBUG
    char residueString[81];
    stringifyResidue(residue, residueString);
    printf ("Destroying residue %s\n", residueString);
#endif
    free(residue->contactResidueList);
    free(residue->resName);
    free(residue->resID);
    if(residue->ext_chainID != NULL)
        free(residue->ext_chainID);
    free(residue);
    return NULL;
}

atom_t *destroyAtomList(atom_t *atomList, int nAtom) {
    //atom_t *next_atom = NULL;
#ifdef DEBUG
    fprintf (stderr, "Destroying atomList\n");
#endif

    for (int n = 0 ; n < nAtom ; n++)
        clearAtomFields(&atomList[n]);
    free(atomList); // Destroy the array

#ifdef DEBUG
    printf("All %d atoms list destroyed\n", nAtom);
#endif
    return NULL;
}
atom_t *destroyAtom(atom_t *atom) {
    clearAtomFields(atom);
    free(atom);
    return NULL;
}
void clearAtomFields(atom_t *atom){
#ifdef DEBUG
    char atomString[81];
    stringifyAtom(atom, atomString);
    fprintf(stderr, "Clearing atom %s\n", atomString);
#endif

    free(atom->resID);
    free(atom->resName);
    free(atom->name);
    if (atom->f_grid != NULL)
        destroyFiboGrid(atom->f_grid);
    if (atom->ext_chainID != NULL)
        free(atom->ext_chainID);
    //free(atom); No need to free atom structure itself, malloc was operate on an array
#ifdef DEBUG
    fprintf(stderr, "Ok\n");
#endif
}
unsigned int atomPairListLen(atomPair_t *atomPairList){
    unsigned int n = 0;
    atomPair_t *currAtomPair = atomPairList;
    while(currAtomPair != NULL) {
        n++;
        currAtomPair = currAtomPair->next;
    }
    return n;
}
unsigned int atomListLen(atom_t *atomList) {
    unsigned int n = 0;
    atom_t *currAtom = atomList;
    while(currAtom != NULL) {
        n++;
        currAtom = currAtom->nextAtomList;
    }
    return n;
}
// Here we add 
atom_t *CreateAtomListFromPdbContainer(pdbCoordinateContainer_t *pdbCoordinateContainer, int *nAtom, atom_map_t *aMap, float probeRadius, int sasaResLvl) {
    #ifdef DEBUG
    fprintf(stderr, "Entering CreateAtomListFromPdbContainer bASA: %s probe radius:%g\n", aMap != NULL ?"true":"false",probeRadius);
    #endif
    double *x;
    double *y;
    double *z;
    char *chainID;
    char **resSeq;
    char **resName;
    char **atomName;
    //atom_t *atomList = NULL;
    *nAtom = pdbContainerToArrays(pdbCoordinateContainer, &x, &y, &z, &chainID, &resSeq, &resName, &atomName);
    atom_t *atomList = readFromArrays(*nAtom, x, y, z, chainID, resSeq, resName, atomName, aMap, probeRadius, sasaResLvl);
    
    freeAtomListCreatorBuffers(x, y, z, chainID, resSeq, resName, atomName, *nAtom);
    #ifdef DEBUG
    fprintf(stderr, "Exiting CreateAtomListFromPdbContainer\n");
    #endif
    return atomList;
}

void freeAtomListCreatorBuffers(double *x, double *y, double *z, char *chainID,\
                                char **resID, char **resName,  char **name, int n) {
#ifdef DEBUG
fprintf(stderr, "freeAtomListCreatorBuffers over even-sized buffers of %d items\n", n); 
#endif
    free(x);
    free(y);
    free(z);
    free(chainID);
    for (int i = 0; i < n ; i++) {
        free(resID[i]);
        free(resName[i]);
        free(name[i]);
    }
    free(resID);
    free(resName);
    free(name);
}

// A simple nbFrame x nbAtoms 2_uple storing individual ASA and frac foreach residues along the frame
sasaFrame_t *createSasaFrame(atom_t *atomList, int nbFrame) {
    sasaFrame_t *sasaFrame = malloc(sizeof(sasaFrame_t));
    sasaFrame->nbFrame = nbFrame;
    sasaFrame->residueList = createResidueList(atomList);
    sasaFrame->nbRes = sasaFrame->residueList->length; // or -1 :p
    
    sasaFrame->sasa2upleArray = malloc(sasaFrame->nbFrame * sizeof(sasa2uple_t *));
    for (int i = 0 ; i < sasaFrame->nbFrame ; i++ )
        sasaFrame->sasa2upleArray[i] = malloc( sasaFrame->nbRes * sizeof(sasa2uple_t) );

    return sasaFrame;
}

sasaFrame_t *destroySasaFrame(sasaFrame_t *sasaFrame) {
    for (int i = 0 ; i < sasaFrame->nbFrame ; i++ )
            free(sasaFrame->sasa2upleArray[i]);
    free(sasaFrame->sasa2upleArray);
    destroyResidueList(sasaFrame->residueList);
    free(sasaFrame);
    return NULL;
}   


//https://numpy.org/doc/stable/user/c-info.how-to-extend.html#example
#ifdef AS_PYTHON_EXTENSION

/* Read multi coordinates atom specs */
atom_t *readFromNumpyArraysFrame(coorFrame_t **coorFrame, PyObject *positionFrame,\
                                PyArrayObject *names, PyArrayObject *resnames,\
                                PyArrayObject *resids, PyArrayObject *segids,\
                                atom_map_t *aMap, float probeRadius, int resolutionLevel) {
    
    // We copy all coordinates frame
    *coorFrame = createFrameFromPyArrayList(positionFrame);
    // We read-in 1st coordinates again to create the atom_t list
    PyObject *firstPositions = PyList_GetItem(positionFrame, 0);
   
    Py_INCREF(firstPositions);
    assert(PyArray_Check(firstPositions));
    atom_t *atomList = readFromNumpyArrays((PyArrayObject*)firstPositions, names, resnames, resids, segids, aMap, probeRadius, resolutionLevel);
    Py_DECREF(firstPositions);
    return atomList;
}

/*
Copy a list of numpy "positions" arrays into a native C structure
*/

coorFrame_t *destroyCoorFrame(coorFrame_t *coorFrame, int optIndex) {
    int nbFrame2Del = optIndex == -1 ? coorFrame->nbFrame : optIndex;
    for (int i = 0; i < nbFrame2Del ; i++) {
        free(coorFrame->coordinates[i]);
    }
    free(coorFrame->coordinates);
    free(coorFrame);
    return NULL;
}
//https://stackoverflow.com/questions/26384129/calling-c-from-python-passing-list-of-numpy-pointers
coorFrame_t *createFrameFromPyArrayList(PyObject *positionArrayList) {
    
    PyArrayObject *positionsBuffer = NULL;
    PyObject *x, *y, *z;
   
    coorFrame_t *coorFrame = malloc(sizeof(coorFrame_t));
    coorFrame->nbFrame     = (int)PyList_Size(positionArrayList);
    coorFrame->nbAtom      = -1;
    coorFrame->coordinates = malloc( coorFrame->nbFrame * sizeof(coordinates_t *) );

    int cLen = 0;
    for (int iFrame = 0 ; iFrame < coorFrame->nbFrame ; iFrame++) {
        positionsBuffer = (PyArrayObject *)PyList_GetItem(positionArrayList, iFrame);
        Py_INCREF(positionsBuffer);
    
        
        npy_intp *shapes = PyArray_SHAPE(positionsBuffer);
        cLen = shapes[0];
        if(coorFrame->nbAtom == -1)
            coorFrame->nbAtom = cLen;
        if(coorFrame->nbAtom != cLen) {
            coorFrame = destroyCoorFrame(coorFrame, iFrame - 1 );
            PySys_WriteStderr("Fatal: Missmatch in coordinates array sizes\n");
            break;
        }

        coorFrame->coordinates[iFrame] = malloc( coorFrame->nbAtom * sizeof(coordinates_t) );
        for (int iAtom = 0 ; iAtom < coorFrame->nbAtom ; iAtom++) {
            //fprintf(stderr, "[%d/%d] %d / %d\n", iFrame, coorFrame->nbFrame, iAtom, coorFrame->nbAtom);
            x = PyArray_GETITEM(positionsBuffer, PyArray_GETPTR2(positionsBuffer, iAtom, 0) );
            y = PyArray_GETITEM(positionsBuffer, PyArray_GETPTR2(positionsBuffer, iAtom, 1) );
            z = PyArray_GETITEM(positionsBuffer, PyArray_GETPTR2(positionsBuffer, iAtom, 2) );
        /*
            fprintf(stderr, "\t[%d] %f %f %f\n",iAtom,\
            (float)PyFloat_AS_DOUBLE(x), (float)PyFloat_AS_DOUBLE(y), (float)PyFloat_AS_DOUBLE(z));
        */
            coorFrame->coordinates[iFrame][iAtom].x = (float)PyFloat_AS_DOUBLE(x);
            coorFrame->coordinates[iFrame][iAtom].y = (float)PyFloat_AS_DOUBLE(y);
            coorFrame->coordinates[iFrame][iAtom].z = (float)PyFloat_AS_DOUBLE(z);
            //fprintf(stderr, "\t Decrefin\n");
            
            Py_DECREF(x);
            Py_DECREF(y);
            Py_DECREF(z);
        }
        //fprintf(stderr, "iFrame %d Done??\n", iFrame);
        Py_DECREF(positionsBuffer);
    }

    return coorFrame;
}

atom_t *readFromNumpyArrays(PyArrayObject *_positions, PyArrayObject *_names,\
                            PyArrayObject *_resnames,  PyArrayObject *_resids, PyArrayObject *_segids,\
                            atom_map_t *aMap, float probeRadius, int resolutionLevel){

    npy_intp *shapes = PyArray_SHAPE(_names);
    int nAtoms = shapes[0];
   
    atom_t *atomList = malloc(nAtoms * sizeof(atom_t));
   
    PyObject *name, *x, *y, *z, *resid, *resname, *segid = NULL;
    char buffer[81];
   
    int resid_buf;
    for (int i = 0 ; i < nAtoms ; i++) {

        name    = PyArray_GETITEM(_names, PyArray_GETPTR1(_names, i) );
        x       = PyArray_GETITEM(_positions, PyArray_GETPTR2(_positions, i, 0) );
        y       = PyArray_GETITEM(_positions, PyArray_GETPTR2(_positions, i, 1) );
        z       = PyArray_GETITEM(_positions, PyArray_GETPTR2(_positions, i, 2) );
        resid   = PyArray_GETITEM(_resids, PyArray_GETPTR1(_resids, i) );
        resname = PyArray_GETITEM(_resnames, PyArray_GETPTR1(_resnames, i) );
        segid   = PyArray_GETITEM(_segids, PyArray_GETPTR1(_segids, i) );

        atomList[i].index = i;
        atomList[i].nextAtomList = NULL;
        atomList[i].belongsTo = NULL;
        atomList[i].inCell = NULL;
        atomList[i].nextResidueAtom = NULL;

        // atom name
        PyObject_ToChar(name,  buffer);
        atomList[i].name = malloc( (strlen(buffer) + 1) * sizeof(char));
        strcpy(atomList[i].name, buffer);
        Py_DECREF(name);

        // atom resid, provided as integer by mdanalysis    
        //if(PyArray_IsIntegerScalar(resid)) 
        sprintf(buffer, "%d", (int)PyLong_AsLong(resid) );
        atomList[i].resID = malloc( (strlen(buffer) + 1) * sizeof(char));
        strcpy(atomList[i].resID, buffer);
        Py_DECREF(resid);

        // atom resname
        PyObject_ToChar(resname,  buffer);
        atomList[i].resName = malloc( (strlen(buffer) + 1) * sizeof(char));
        strcpy(atomList[i].resName, buffer);
        Py_DECREF(resname);
        
        // atom segid, a string in MDanalysis
        PyObject_ToChar(segid,  buffer);
        atomList[i].ext_chainID = malloc( (strlen(buffer) + 1) * sizeof(char));
        strcpy(atomList[i].ext_chainID, buffer);
        Py_DECREF(segid);
            // Setting native chainID to space char
        atomList[i].chainID = ' ';

        // coordinates
        atomList[i].x = (float)PyFloat_AS_DOUBLE(x);
        atomList[i].y = (float)PyFloat_AS_DOUBLE(y);
        atomList[i].z = (float)PyFloat_AS_DOUBLE(z);
        Py_DECREF(x);
        Py_DECREF(y);
        Py_DECREF(z);
        /*
        if(i%1000 == 0) {
            fprintf(stderr, ">%s<\n", atomList[i].name);
            fprintf(stderr, ">%s<\n", atomList[i].resID);
            fprintf(stderr, ">%s<\n", atomList[i].resName);
            fprintf(stderr, ">%s<\n", atomList[i].ext_chainID);
            fprintf(stderr, ">%f %f %f<\n",  atomList[i].x,\
            atomList[i].y,  atomList[i].z);
        }
        */
        if (i > 0)
            atomList[i - 1].nextAtomList = &atomList[i];
        
        atomList[i]._radiusASA = aMap != NULL\
                                ? getRadius(aMap, atomList[i].name, atomList[i].resName)\
                                : VDW_DEFAULT; 
        atomList[i]._radiusASA += probeRadius;

        atomList[i].f_grid = aMap != NULL\
            ? computeFiboGrid(atomList[i].x, atomList[i].y, atomList[i].z, atomList[i]._radiusASA, resolutionLevel)\
            : NULL;
        
    }
    return atomList;
}
#endif

atom_t *createBareboneAtom(int n, double x, double y, double z, char chainID, char *resID, \
                           char *resName, char *name) {
        atom_t *atom = malloc( 1 * sizeof(atom_t));
        atom->index = n;
        
        atom->resID = malloc( (strlen(resID) + 1) * sizeof(char));
        strcpy(atom->resID, resID);

        atom->resName = malloc( (strlen(resName) + 1) * sizeof(char));
        strcpy(atom->resName, resName);

        atom->name = malloc( (strlen(name) + 1) * sizeof(char));
        strcpy(atom->name, name);

        atom->nextAtomList = NULL;
        atom->belongsTo = NULL;
        atom->inCell = NULL;
        atom->nextResidueAtom = NULL;
        atom->x = x;
        atom->y = y;
        atom->z = z;
        atom->chainID     = chainID;
        atom->f_grid = NULL;
        atom->ext_chainID = NULL;
        atom->_VDWradius = 0.0;
        atom->_radiusASA = 0.0;

    return atom;
}
// MEMORY ALLOCATION OF atom LIST
atom_t *readFromArrays(int nAtoms, double *x, double *y, double *z, char *chainID, \
                       char **resID, char **resName, char **name, atom_map_t *aMap,\
                       float probeRadius, int resolutionLevel) { // Resume here
    #ifdef DEBUG
        char DBG_buffer[200];
        sprintf(DBG_buffer, "readFromArray: Running readFromArrays over %d atoms\n", nAtoms);
        printOnContextStderr(DBG_buffer);
        for (int i = 0 ; i < nAtoms ; i++) {
            sprintf(DBG_buffer, "[%g, %g, %g] %c, %s, %s %s\n", x[i], y[i], z[i], chainID[i], resID[i], resName[i], name[i]);
            printOnContextStderr(DBG_buffer);
        }
    #endif
   
    atom_t *atomList = malloc(nAtoms * sizeof(atom_t));
    for (int n = 0 ; n < nAtoms ; n++) {
        atomList[n].index = n;
        atomList[n].nextAtomList = NULL;
        atomList[n].belongsTo = NULL;
        atomList[n].inCell = NULL;
        atomList[n].nextResidueAtom = NULL;

        atomList[n].resID = malloc( (strlen(resID[n]) + 1) * sizeof(char));
        strcpy(atomList[n].resID, resID[n]);

        atomList[n].resName = malloc( (strlen(resName[n]) + 1) * sizeof(char));
        strcpy(atomList[n].resName, resName[n]);

        atomList[n].name = malloc( (strlen(name[n]) + 1) * sizeof(char));
        strcpy(atomList[n].name, name[n]);

        atomList[n].x = x[n];
        atomList[n].y = y[n];
        atomList[n].z = z[n];
        atomList[n].chainID     = chainID[n];
        atomList[n].ext_chainID = NULL;

        if (n > 0)
            atomList[n - 1].nextAtomList = &atomList[n];
       
        atomList[n]._VDWradius = aMap != NULL ? getRadius(aMap, atomList[n].name, atomList[n].resName)\
                                      : VDW_DEFAULT; 
        atomList[n]._radiusASA = atomList[n]._VDWradius; 
        atomList[n]._radiusASA += probeRadius;
        #ifdef DEBUG
            sprintf(DBG_buffer, "Assiging to atom object[%g, %g, %g] %c, %s, %s %s %g\n",\
                atomList[n].x, atomList[n].y, atomList[n].z, atomList[n].chainID, atomList[n].resID,\
                atomList[n].resName, atomList[n].name, atomList[n]._VDWradius, atomList[n]._radiusASA);
            printOnContextStderr(DBG_buffer);
        #endif
        
        if (resolutionLevel > -1) {           
            atomList[n].f_grid = computeFiboGrid(atomList[n].x, atomList[n].y, atomList[n].z, atomList[n]._radiusASA, resolutionLevel);
            #ifdef DEBUG
            printOnContextStderr("f_grid build succesfull\n");
            #endif
        } else {
            atomList[n].f_grid = NULL;
        }



    }
    
    #ifdef DEBUG
        sprintf(DBG_buffer, "Read from array done for %d\nTrying to show atom_t *atomList content\n",\
                            nAtoms);
        printOnContextStderr(DBG_buffer);
        
        //printAtomList(atomList, stderr);
        char *atomListCharList = stringifyAtomList(atomList);
        printOnContextStderr(atomListCharList);    
        free(atomListCharList);
        printOnContextStderr("\n--Exiting readFromArrays--\n");
    #endif

    return atomList;
}

atom_t *legacy_readCoordinates(char *fname, int *_nAtom) {
    atom_t *head = NULL;
    atom_t *old = NULL;
    head = malloc(sizeof(atom_t));
    atom_t *root = head;
    FILE * fp;
    size_t len = 0;
    size_t read;
    char * line = NULL;

    fp = fopen(fname, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    char *p;
    char *end;

    int nAtom = 0;
    while ((int)(read = getline(&line, &len, fp)) != -1) {
        p = line;
        double buf[3];
        int i = 0;
        for (double f = strtod(p, &end); p != end; f = strtod(p, &end)) {
            buf[i] = f;
            p = end;
            //printf("%f\n", f);
            i++;
        }

        head->x = buf[0];
        head->y = buf[1];
        head->z = buf[2];
        head->index = 0;
        head->nextAtomList = malloc(sizeof(atom_t));
        head = head->nextAtomList;
        head->nextAtomList = NULL;
        nAtom++;
    }
    fclose(fp);
    if (line)
        free(line);
    atom_t *atomArray = malloc (nAtom * sizeof(atom_t));
    head = root;
    int i = 0;
    while (head != NULL) {
        if (i < nAtom) {
            atomArray[i].x = head->x;
            atomArray[i].y = head->y;
            atomArray[i].z = head->z;
        }
        old = head;
        head = head->nextAtomList;
        free(old);
        i++;
    }
    free(head);

    *_nAtom = nAtom;
    return atomArray;
}
int appendFiboGridToPdbContainer(pdbCoordinateContainer_t *cloudPdbContainer,\
                                atom_t *atomList, int nbAtom, char segID, bool hiddenBuried)
{
    atom_t *currAtom = atomList;
    if (currAtom == NULL) {
        fprintf(stderr, "appendFiboGridToPdbContainer: Unexpected NULL atom pointer at atom list head\n");
        return-1;
    }
    fibo_grid_t *fibo_grid = currAtom->f_grid;
    if (fibo_grid == NULL) {
        fprintf(stderr, "appendFiboGridToPdbContainer: Unexpected NULL fibonacci grid pointer at atom list head\n");
        return-1;
    }
    /* --  Setting coordinates total arrays and buffers -- */
    int fiboGridSize = currAtom->f_grid->n_spots;
    // Estimating realloc chunk size as 100 * total number of points on a single fibonacci grid
    //int coorChunk = 100 * (size_t)fiboGridSize;
    int chunk_sz = 100;
    int nb_chunks = 1;
    int capacity = nb_chunks * chunk_sz;
    double *x_total = malloc(capacity * (sizeof(double)) );
    double *y_total = malloc(capacity * (sizeof(double)) );
    double *z_total = malloc(capacity * (sizeof(double)) );
   

    double *x_buffer = malloc(sizeof(double) * fiboGridSize);
    double *y_buffer = malloc(sizeof(double) * fiboGridSize);
    double *z_buffer = malloc(sizeof(double) * fiboGridSize);
    /* -- -- */

    // Generate remaining arrays for pdb container
    char baseResNameV[4] = "VOX" ;
    char baseResNameS[4] = "SOX" ;

   
    char **resName_total  = malloc(capacity * sizeof(char*));
    char **resID_total    = malloc(capacity * sizeof(char*));
    char **atomName_total = malloc(capacity * sizeof(char*));
    for(int  i = 0; i < capacity ; i++) {
        resName_total[i]  = malloc(sizeof(char) * 4);
        atomName_total[i] = malloc(sizeof(char) * 5);
        resID_total[i]    = malloc(sizeof(char) * 5);
    }
    char *segID_total    = malloc(capacity * sizeof(char));


    int curPgridCount;
    int totPgridCount = 0;
    int prev_capacity;
    while(currAtom != NULL) {
        fibo_grid = currAtom->f_grid;
        curPgridCount = generateGridPointCartesian(fibo_grid, x_buffer, y_buffer, z_buffer, hiddenBuried);
       // fprintf(stderr, "Testing list access %f %f %f [%d/%d]\n",
       //         x_buffer[0], y_buffer[0], z_buffer[0], 0, curPgridCount);

        if (totPgridCount + curPgridCount > capacity) {          
           /*
            fprintf(stderr, "Reallocating fibo grid points coordinates from %d to %d\n",\
                capacity, (nb_chunks+1) * chunk_sz);
            */
            // realloc
            nb_chunks++;
            capacity = nb_chunks * chunk_sz;
            x_total = realloc(x_total, capacity * sizeof(double));
            y_total = realloc(y_total, capacity * sizeof(double));
            z_total = realloc(z_total, capacity * sizeof(double));
            resName_total  = realloc(resName_total, capacity * sizeof(char*));
            resID_total    = realloc(resID_total,   capacity * sizeof(char*));
            atomName_total = realloc(atomName_total, capacity * sizeof(char*));
            
            prev_capacity  = (nb_chunks-1) * chunk_sz;
            for(int  i = prev_capacity; i < capacity ; i++) {
                resName_total[i]  = malloc(sizeof(char) * 4);
                atomName_total[i] = malloc(sizeof(char) * 5);
                resID_total[i]    = malloc(sizeof(char) * 5);
            }
            segID_total    = realloc(segID_total, capacity * sizeof(char));
        }
        // Copy
        for (int i = 0 ; i < curPgridCount ; i++) {
            x_total[totPgridCount + i] = x_buffer[i];
            y_total[totPgridCount + i] = y_buffer[i];
            z_total[totPgridCount + i] = z_buffer[i];
            strcpy(resName_total[totPgridCount + i], \
                hiddenBuried ? baseResNameS : baseResNameV);
            strcpy(resID_total[totPgridCount + i], currAtom->resID);
            strcpy(atomName_total[totPgridCount + i], currAtom->name);
            segID_total[totPgridCount + i] = segID;
        }

        totPgridCount += curPgridCount;
        currAtom = currAtom->nextAtomList;
    }

    free(x_buffer);
    free(y_buffer);
    free(z_buffer);

    appendArraysToPdbContainer(cloudPdbContainer, totPgridCount,\
        x_total, y_total, z_total, segID_total, resID_total, resName_total, atomName_total);

    for(int  i = 0; i < capacity ; i++) {
        free(resName_total[i]);
        free(atomName_total[i]);
        free(resID_total[i]);
    }
    free(resName_total);
    free(atomName_total);
    free(resID_total);
    free(segID_total);
    free(x_total);
    free(y_total);
    free(z_total);

    return totPgridCount;
}
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
    char atomString[81];
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
    

    char atomString[81];
    #ifdef DEBUG
        char dbgBuffer[81];
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
    sprintf(atomString, "%s %s %s %c %f %f %f", atom->name, atom->resName, atom->resID, atom->chainID, atom->x, atom->y, atom->z);    
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

// Replace the x,yz values of the 2nd list by those of the 1st list
// Compute fibo grid if needed
// Returns true if lists were of even size  -> FIBO GRID RADIAUS VALUE IS WRONG
bool applyCoordinates(atom_t *atomListFrom, atom_t *atomListTo) {
    
    while (atomListFrom != NULL && atomListTo != NULL ) {
        atomListTo->x = atomListFrom->x;
        atomListTo->y = atomListFrom->y;
        atomListTo->z = atomListFrom->z;
        if(atomListFrom->f_grid != NULL) 
           atomListTo->f_grid = computeFiboGrid(atomListTo->x, atomListTo->y, atomListTo->z, atomListTo->_radiusASA);
        
        atomListFrom = atomListFrom->nextAtomList;
        atomListTo   = atomListTo->nextAtomList;
    }
    
    return atomListTo == NULL && atomListFrom == NULL;
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

    while(curr_atom != NULL) {
        if (strcmp(curr_atom->resID, residue_head->resID) == 0
            && curr_atom->chainID == residue_head->chainID) {
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
    free(residue);
    return NULL;
}

atom_t *destroyAtomList(atom_t *atomList, int nAtom) {
    //atom_t *next_atom = NULL;
#ifdef DEBUG
    printf ("Destroying atomList\n");
#endif

    for (int n = 0 ; n < nAtom ; n++)
        destroyAtom(&atomList[n]);
    free(atomList); // Destroy the array

#ifdef DEBUG
    printf("All %d atoms list destroyed\n", nAtom);
#endif
    return NULL;
}

atom_t *destroyAtom(atom_t *atom){
#ifdef DEBUG
    char atomString[81];
    stringifyAtom(atom, atomString);
    printf ("Destroying atom %s\n", atomString);
#endif

    free(atom->resID);
    free(atom->resName);
    free(atom->name);
    if (atom->f_grid != NULL)
        destroyFiboGrid(atom->f_grid);
    //free(atom); No need to free atom structure itself, malloc was operate on an array
#ifdef DEBUG
    printf("Ok\n");
#endif
return NULL;
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

atom_t *CreateAtomListFromPdbContainer(pdbCoordinateContainer_t *pdbCoordinateContainer, int *nAtom, atom_map_t *aMap, float probeRadius) {
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
    atom_t *atomList = readFromArrays(*nAtom, x, y, z, chainID, resSeq, resName, atomName, aMap, probeRadius);
    
    freeAtomListCreatorBuffers(x, y, z, chainID, resSeq, resName, atomName, *nAtom);
    #ifdef DEBUG
    fprintf(stderr, "Exiting CreateAtomListFromPdbContainer\n");
    #endif
    return atomList;
}

void freeAtomListCreatorBuffers(double *x, double *y, double *z, char *chainID, char **resID, char **resName,  char **name, int n) {
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

//https://numpy.org/doc/stable/user/c-info.how-to-extend.html#example

atom_t *readFromNumpyArrays(PyArrayObject *_positions, PyArrayObject *_names,\
                            PyArrayObject *_resnames,  PyArrayObject *_resids, PyArrayObject *_segids,\
                            atom_map_t *aMap, float probeRadius){

    int nAtoms = (int) PyArray_SIZE(_names);
   
    atom_t *atomList = malloc(nAtoms * sizeof(atom_t));
    PySys_WriteStderr(">>%d<<\n", nAtoms);
    PyObject *name, *x, *y, *z, *resid, *resname, *segid = NULL;
    char buffer[81];
   
    int resid_buf;
    for (int i = 0 ; i < nAtoms ; i++) {

        name    = PyArray_GETITEM(_names, PyArray_GETPTR1(_names, i) );
     //   x       = PyArray_GETITEM(_positions, PyArray_GETPTR2(_positions, i, 0) );
     //   y       = PyArray_GETITEM(_positions, PyArray_GETPTR2(_positions, i, 1) );
     //   z       = PyArray_GETITEM(_positions, PyArray_GETPTR2(_positions, i, 2) );
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

        // atom resid, provided as integer by mdanalysis
        // int my_val = (int) PyInt_AsLong(pyObj_val);
        sprintf(buffer, "%l", PyLong_AsLong(resid) );
        atomList[i].resID = malloc( (strlen(buffer) + 1) * sizeof(char));
        strcpy(atomList[i].resID, buffer);

        if(i%1000 == 0) {
            fprintf(stderr, "%s\n", atomList[i].name);
            fprintf(stderr, "%s\n", atomList[i].resID);
        }
       // PyObject_ToChar(resname,  buffer);
       /*
        PyObject_ToChar(resid,  buffer);
        if(i%1000 == 0) 
            fprintf(stderr, "->%s\n", buffer);

       // atomList[i].resID = malloc( (strlen(buffer) + 1) * sizeof(char));
       // strcpy(atomList[i].resID, buffer);
        */
      /*  if(i%1000 == 0)
            fprintf(stderr, "%s||%s\n", atomList[i].name, atomList[i].resID\
            );*/

/*
        PyObject_ToChar(resname,  buffer);
        atomList[i].resName = malloc( (strlen(buffer) + 1) * sizeof(char));
        strcpy(atomList[i].resName, buffer);

        PyObject_ToChar(segid,  buffer);
        atomList[i].chainID = buffer[0];

        atomList[i].x = (float)PyFloat_AS_DOUBLE(x);
        atomList[i].y = (float)PyFloat_AS_DOUBLE(y);
        atomList[i].z = (float)PyFloat_AS_DOUBLE(z);
       
        stringifyAtom(&(atomList[i]), buffer);
        PySys_WriteStdout("%s\n", buffer);

        */
        /*
        PySys_WriteStdout("%s\n", buffer);
        PyObject_ToChar(name,  name_buffer);
        */
        /*
        PySys_WriteStderr("### %s name %f %f %f\n",\
            name_buffer, *x, *y, *z);//names[i]);
        */
       
        Py_DECREF(name);
        Py_DECREF(resid);

       /* Py_DECREF(x);
        Py_DECREF(y);
        Py_DECREF(z);
  
        Py_DECREF(resname);
        Py_DECREF(segid);
        */
    }
    return NULL;//atomList;
}

/*    atom_t *atomList = malloc(nAtoms * sizeof(atom_t));
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
        atomList[n].chainID = chainID[n];

        if (n > 0)
            atomList[n - 1].nextAtomList = &atomList[n];

        atomList[n]._radiusASA = aMap != NULL ? getRadius(aMap, atomList[n].name, atomList[n].resName)\
                                           : VDW_DEFAULT; 
        atomList[n]._radiusASA += probeRadius;
        
        #ifdef DEBUG
            sprintf(DBG_buffer, "Assiging to atom object[%g, %g, %g] %c, %s, %s %s %g\n",\
                atomList[n].x, atomList[n].y, atomList[n].z, atomList[n].chainID, atomList[n].resID,\
                atomList[n].resName, atomList[n].name, atomList[n]._radiusASA);
            printOnContextStderr(DBG_buffer);
        #endif
        
        if (aMap != NULL) {           
            atomList[n].f_grid = computeFiboGrid(atomList[n].x, atomList[n].y, atomList[n].z, atomList[n]._radiusASA);
            #ifdef DEBUG
            printOnContextStderr("f_grid build succesfull\n");
            #endif
        } else {
            atomList[n].f_grid = NULL;
        }



    }
    
    */

//}

// MEMORY ALLOCATION OF atom LIST
atom_t *readFromArrays(int nAtoms, double *x, double *y, double *z, char *chainID, char **resID, char **resName, char **name, atom_map_t *aMap, float probeRadius) { // Resume here
    #ifdef DEBUG
        char DBG_buffer[81];
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
        atomList[n].chainID = chainID[n];

        if (n > 0)
            atomList[n - 1].nextAtomList = &atomList[n];

        atomList[n]._radiusASA = aMap != NULL ? getRadius(aMap, atomList[n].name, atomList[n].resName)\
                                           : VDW_DEFAULT; 
        atomList[n]._radiusASA += probeRadius;
        
        #ifdef DEBUG
            sprintf(DBG_buffer, "Assiging to atom object[%g, %g, %g] %c, %s, %s %s %g\n",\
                atomList[n].x, atomList[n].y, atomList[n].z, atomList[n].chainID, atomList[n].resID,\
                atomList[n].resName, atomList[n].name, atomList[n]._radiusASA);
            printOnContextStderr(DBG_buffer);
        #endif
        
        if (aMap != NULL) {           
            atomList[n].f_grid = computeFiboGrid(atomList[n].x, atomList[n].y, atomList[n].z, atomList[n]._radiusASA);
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

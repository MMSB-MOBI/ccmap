//File: pdb_coordinates.h
#ifndef PDB_COORDINATES_H
#define PDB_COORDINATES_H

typedef struct atomRecord {
    char recordName[7];
    int serial;
    char name[5];
    char altLoc;
    char resName[4];
    char chainID;
    char resSeq[5];
    char iCode;
    double x;
    double y;
    double z;
    double occupancy;
    double tempFactor;
    char element[3];
    char charge[3];
} atomRecord_t;

typedef struct pdbCoordinateContainer {
    struct atomRecord *atomRecordArray;
    int atomCount;
} pdbCoordinateContainer_t;
/*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/

void transformPdbCoordinateContainer(pdbCoordinateContainer_t *pdbCoordinateContainer, float *euler, float *translation);

void pdbContainerToFile(pdbCoordinateContainer_t *pdbCoordinateContainer, char *fname, char *mode);
pdbCoordinateContainer_t *pdbFileToContainer(char *fileName);
pdbCoordinateContainer_t *destroyPdbCoordinateContainer(pdbCoordinateContainer_t *pdbCoordinateContainer);
void createAtomRecord(char *recordString, atomRecord_t *newAtom);
void stringifyAtomRecord(atomRecord_t *atomRecord, char *atomRecordString);
int pdbContainerToArrays(pdbCoordinateContainer_t *pdbCoordinateContainer, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name);
char *pdbContainerToString(pdbCoordinateContainer_t *pdbCoordinateContainer);


// Unused kept just in case
int legacy_readPdbFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name);
int legacy_readFile(char *fname, double **x, double **y, double **z, char **chainID, char ***resID, char ***resName,  char ***name);
#endif

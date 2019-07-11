typedef enum { false, true } bool;

typedef struct atomPair {
    struct atom *a;
    struct atom *b;
    double dist;
    struct atomPair *next;
} atomPair_t;

typedef struct meshContainer {
    struct mesh *mesh;
    struct cell **filledCells;
    int nFilled;
    float step; // Grid step
} meshContainer_t;

typedef struct residue {
    struct atom *elements;
    int nAtoms;
    int index;
    char *resID;
    char chainID;
    char *resName;
    struct residue *nextResidueList;
    struct residue *prevResidueList;
    struct residue **contactResidueList; // list of residue ptr dynmically incremented
    int nContacts;
   // struct residue
} residue_t;

typedef struct atom {
    struct residue *belongsTo;
    struct cell *inCell;
    char chainID;
    char *resID;
    float x;
    float y;
    float z;
    char *resName;
    char *name;
    struct atom *nextAtomList;
    struct atom *nextCellAtom;
    struct atom *nextResidueAtom;
} atom_t;

typedef struct cell {
    struct atom *members;
    struct atom *head;
    int memberCount;

    struct atom *iMembers;
    struct atom *iHead;
    int iMemberCount;

    struct atom *jMembers;
    struct atom *jHead;
    int jMemberCount;


    int i;
    int j;
    int k;
    int n;
    struct neighbours *cell;
    int neighbourCount;
} cell_t;

typedef struct mesh {
    int iMax;
    int jMax;
    int kMax;
    int n;
    cell_t ***grid;
} mesh_t;

//Utility functions
int concatenate(char **dest, char *src);
int popChar(char **dest, char c);
//DEBUG
void dumpBuffer(char *buffer, int bufSize);

// this is the inputs processing function
atom_t *readFromArrays(int nAtoms, double *x, double *y, double *z, char *chainID, char **resSeq, char **resName, char **atomName);
// this is the main computing function
char *residueContactMap(atom_t * atomList, int nAtom, double ctc_dist);
char *residueContactMap_DUAL(atom_t *iAtomList, int iNatom, atom_t *jAtomList, int jNatom, double ctc_dist);

// Display structure content
void printResidueList(residue_t *residueList);
void printResidue(residue_t *residue);
void printAtomList(atom_t *atomList);
void stringifyAtom(atom_t *atom, char *atomString);
void stringifyResidue(residue_t *residue, char *residueString);

//Results display
void dumpMeshContent(meshContainer_t *meshContainer);
void dumpCellContent(cell_t *cell);
void printContactList(residue_t *residueList);
char *jsonifyContactList(residue_t *residueList);
void jsonifyResidue(residue_t *residue, char *jsonString);
void atomListInContact(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double ctc_dist, int iAtomStatus[], int jAtomStatus[]);



residue_t *createResidue(atom_t *atom, int n);
residue_t *createResidueList(atom_t * atomList);
void fuseResidueLists(residue_t *iResidueList, residue_t *jResidueList);
residue_t *destroyResidueList(residue_t *residueList);
residue_t *destroyResidue(residue_t *residue);
atom_t *destroyAtomList(atom_t *atomList, int nAtom);
atom_t *destroyAtom(atom_t *atom);

int updateContactList(atom_t *iAtom, atom_t *jAtom);
int updateContactList_DUAL(atom_t *iAtom, atom_t *jAtom);

mesh_t *createMesh(int iDim, int jDim, int kDim);
meshContainer_t *destroyMeshContainer(meshContainer_t *container);
mesh_t *destroyMesh(mesh_t *i_mesh);

void printMesh(mesh_t *mesh);
void enumerate(meshContainer_t *meshContainer, double ctc_dist, int *nPairs, bool dualBool);
void pairwiseCellEnumerate(cell_t *refCell, cell_t *targetCell, double ctc_dist, int *nContacts, int *nDist);
void pairwiseCellEnumerate_DUAL(cell_t *refCell, cell_t *targetCell, double ctc_dist, int *nContacts, int *nDist);
void getBoundariesCartesian(atom_t * atomList, int nAtom, atom_t *minCoor, atom_t *maxCoor);
void cartesianToMesh(atom_t *atom, int *i, int *j, int *k, float step, atom_t minCoor);
meshContainer_t *createMeshContainer(atom_t *iAtomList, int iAtom, atom_t *jAtomList, int jAtom, double step);
//void mesh(atom_t * atomList, int nAtom, double step);
cell_t ** vectorizeMesh(mesh_t *mesh);
void meshDummy(int a, int b, int c);
double distance(atom_t *iAtom, atom_t *jAtom);

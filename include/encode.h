
typedef struct ChainedInt {
    int index;
    struct ChainedInt* nextInt;
} chainedInt_t;

int chainLen(residue_t *ResidueList);
int contactIndex(int index1, int index2, int max2);
int *encodeContactMap(residue_t *ResidueList, int lenLigList, int lenRecList, int *finalLen);
chainedInt_t *createChained_Int(int index);
void printTable(int *ContactList, int len);
void printChain(chainedInt_t *ContactList);
int *copyTable(int *table, int lenTable );

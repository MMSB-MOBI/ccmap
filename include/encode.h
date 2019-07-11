
typedef struct ChainedInt {
    int index;
    struct ChainedInt* nextInt;
} chainedInt_t;

int chainLen(residue_t *ResidueList);
int contactIndex(int index1, int index2, int max2);
chainedInt_t *encodeContactMap(residue_t *ResidueList, int lenLigList);
chainedInt_t *createChained_Int(int index);
void printChain(chainedInt_t *ContactList);

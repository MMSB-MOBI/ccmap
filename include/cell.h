//File: cell.h
#ifndef CELL_H
#define CELL_H

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
#endif

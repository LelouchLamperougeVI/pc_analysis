#ifndef _GET_RANKS
#define _GET_RANKS

struct ranksData {
        int m, idx;
        int *numTasks;
        double *A, *ranks;
};

void *getRanks(void *arg);

#endif

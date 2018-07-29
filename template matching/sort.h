#ifndef _SORT
#define _SORT

struct sortStruct {
        double *value;
        double idx;
};

struct sortData {
        int m, t_idx;
        int *numTasks;
        double *A, *sorted, *idx;
};

void *quicksort(void *arg);

#endif

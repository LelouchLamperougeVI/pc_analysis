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

int cmpfunc(void const *a, void const *b);

void quicksort(struct sortData *data);

#endif

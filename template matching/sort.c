#include "sort.h"

int cmpfunc(void const *a, void const *b){
        struct sortStruct *a1 = (struct sortStruct *) a;
        struct sortStruct *b1 = (struct sortStruct *) b;
        if((double)*(*a1).value < (double)*(*b1).value)
                return -1;
        else if((double)*(*a1).value > (double)*(*b1).value)
                return 1;
        else
                return 0;
}

void quicksort(struct sortData *data){
        int m = data->m;
        int t_idx = data->t_idx;
        int *numTasks = data->numTasks;
        double *idx = data->idx;
        double *A = data->A;
        double *sorted = data->sorted;
        int i, j, block_i;

        struct sortStruct sort[m];

        block_i = 0;
        for(i = 0; i < t_idx; i++)
                block_i += numTasks[i];

        for(i = 0; i < numTasks[t_idx]; i++) {
                for(j = 0; j < m; j++) {
                        sort[j].value = &A[(block_i + i)*m + j];
                        sort[j].idx = j + 1;
                }
                qsort(sort, m, sizeof(struct sortStruct), cmpfunc);
                for(j = 0; j < m; j++) {
                        sorted[(block_i + i)*m + j] = *sort[j].value;
                        idx[(block_i + i)*m + j] = sort[j].idx;
                }
        }
}

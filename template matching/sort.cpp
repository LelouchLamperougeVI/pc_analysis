#include "sort.h"
#include "stdlib.h"
#include <string>
#include <algorithm>

using namespace std;

bool cmpfunc(struct sortStruct a, struct sortStruct b){
        return *a.value < *b.value;
}

void *quicksort(void *arg){
        struct sortData *data = (struct sortData *) arg;
        int m = data->m;
        int t_idx = data->t_idx;
        int *numTasks = data->numTasks;
        double *idx = data->idx;
        double *A = data->A;
        double *sorted = data->sorted;
        int i, j, block_i;

        struct sortStruct sort_arr[m];

        block_i = 0;
        for(i = 0; i < t_idx; i++)
                block_i += numTasks[i];

        for(i = 0; i < numTasks[t_idx]; i++) {
                for(j = 0; j < m; j++) {
                        sort_arr[j].value = &A[(block_i + i)*m + j];
                        sort_arr[j].idx = j + 1;
                }
                sort(sort_arr, sort_arr+m, cmpfunc);
                for(j = 0; j < m; j++) {
                        sorted[(block_i + i)*m + j] = *sort_arr[j].value;
                        idx[(block_i + i)*m + j] = sort_arr[j].idx;
                }
        }
}

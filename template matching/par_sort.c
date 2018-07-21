/** A parallel quicksort implementation
   [sorted, idx] = par_sort(unsorted_matrix)
 */

#include "mex.h"
#include "matrix.h"
#include "pthread.h"
#include "qsort.h" //an inline version of qsort() found here: http://www.corpit.ru/mjt/qsort/qsort.h

#define NUMTHREADS 7

struct sortStruct {
        double *value;
        double idx;
};

struct threadData {
        int m, t_idx;
        int *numTasks;
        double *A, *sorted, *idx;
};

static inline int cmpfunc(void const *a, void const *b){
        struct sortStruct *a1 = (struct sortStruct *) a;
        struct sortStruct *b1 = (struct sortStruct *) b;
        if((double)*(*a1).value < (double)*(*b1).value)
                // return -1;
                return 1;
        // else if((double)*(*a1).value > (double)*(*b1).value)
        //         return 1;
        else
                return 0;
}

void quicksort(struct threadData *data){
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
                // qsort(sort, m, sizeof(struct sortStruct), cmpfunc);
                QSORT(struct sortStruct*, sort, m, cmpfunc);
                for(j = 0; j < m; j++) {
                        sorted[(block_i + i)*m + j] = *sort[j].value;
                        idx[(block_i + i)*m + j] = sort[j].idx;
                }
        }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        if(nrhs > 1)
                mexErrMsgTxt("Too many inputs");
        if(nrhs < 1)
                mexErrMsgTxt("No input given");
        if(nlhs > 2)
                mexErrMsgTxt("Too many outputs assigned");

        double *A, *sorted, *idx;
        int i, m, n;
        int tasksPerThread[NUMTHREADS] = { 0 };

        A = mxGetPr(prhs[0]);
        m = mxGetM(prhs[0]);
        n = mxGetN(prhs[0]);

        plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
        sorted = mxGetPr(plhs[0]);
        plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
        idx = mxGetPr(plhs[1]);

        struct threadData data[NUMTHREADS];
        pthread_t thread[NUMTHREADS];

        for(i = 0; i < n; i++)
                tasksPerThread[i % NUMTHREADS]++;

        for(i = 0; i < NUMTHREADS; i++) {
                data[i].A = A;
                data[i].sorted = sorted;
                data[i].m = m;
                data[i].numTasks = tasksPerThread;
                data[i].idx = idx;
                data[i].t_idx = i;
        }

        for(i = 0; i < NUMTHREADS; i++)
                pthread_create(&thread[i], NULL, quicksort, &data[i]);

        for(i = 0; i < NUMTHREADS; i++)
                pthread_join(thread[i], NULL);

        return;
}

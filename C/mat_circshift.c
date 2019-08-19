/*************************************************************************
 * MATLAB MEX ROUTINE mat_circshift.c
 *
 * Shift row elements of each column by index specified in array
 *
 * Direct complaints to HaoRan Chang
 ************************************************************************/

#include "string.h"
#include "mex.h"
#include "pthread.h"

#define NUMTHREADS 7

struct threadData{
    int idx, m;
    int *numTasks;
    double *A, *shifted, *index;
};

void shiftCol(double *A, double *shifted, int col, int m, int shift){
    int i;
    for(i=0; i<m-shift; i++){
        shifted[col*m + i] = A[col*m + i + shift];
    }
    for(i=m-shift; i<m; i++){
        shifted[col*m + i] = A[col*m + i - m + shift];
    }
}

void *shifter(void *arg){
    struct threadData *data = (struct threadData *) arg;
    int m = data->m;
    double *A = data->A;
    double *index = data->index;
    double *shifted = data->shifted;
    int *numTasks = data->numTasks;
    int idx = data->idx;
    int i;

    int block_i = 0;
    for(i = 0; i < idx; i++)
            block_i += numTasks[i];

    for(i = 0; i < numTasks[idx]; i++){
        shiftCol(A, shifted, block_i+i, m, (int) index[block_i+i]);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if(nrhs != 2){
        if(nrhs < 2){
            mexErrMsgTxt("Two inputs are required");
        }else if(nrhs > 2){
            mexErrMsgTxt("Too many input arguments");
        }
    }

    double *A, *index, *shifted;
    int i, n, m;
    int tasksPerThread[NUMTHREADS] = { 0 };

    A = mxGetPr(prhs[0]);
    index = mxGetPr(prhs[1]);
    n = mxGetN(prhs[1]);
    m = mxGetM(prhs[0]);

    if(mxGetM(prhs[1])>1){
        mexErrMsgTxt("Index must be a columns vector");
    }
    if(n != mxGetN(prhs[0])){
        mexErrMsgTxt("Index and data must have the same number of columns");
    }

    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    shifted = mxGetPr(plhs[0]);

    struct threadData data[NUMTHREADS];
    pthread_t thread[NUMTHREADS];

    for(i = 0; i < n; i++)
            tasksPerThread[i % NUMTHREADS]++;

    for(i=0; i<NUMTHREADS; i++){
        data[i].A = A;
        data[i].shifted = shifted;
        data[i].index = index;
        data[i].m = m;
        data[i].numTasks = tasksPerThread;
        data[i].idx = i;
    }

    for(i=0; i<NUMTHREADS; i++){
        pthread_create(&thread[i], NULL, shifter, &data[i]);
    }

    for(i=0; i<NUMTHREADS; i++){
        pthread_join(thread[i], NULL);
    }

    return;
}

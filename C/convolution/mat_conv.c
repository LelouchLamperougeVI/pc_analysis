/*************************************************************************
 * MATLAB MEX ROUTINE mat_conv.c
 *
 * Simple 1D convolution along column matrices with multithreading
 *
 * Direct complaints to HaoRan Chang
 ************************************************************************/

#include "string.h"
#include "mex.h"
#include "pthread.h"
#include "convolve.h"

#define NUMTHREADS 4

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if(nrhs != 2){
        if(nrhs < 2){
            mexErrMsgTxt("Two inputs are required");
        }else if(nrhs > 2){
            mexErrMsgTxt("Too many input arguments");
        }
    }

    double *A, *kernel, *conv;
    int n, a_size, k_size, i;

    A = mxGetPr(prhs[0]);
    kernel = mxGetPr(prhs[1]);
    n = mxGetN(prhs[0]);
    a_size = mxGetM(prhs[0]);
    k_size = mxGetM(prhs[1]);

    // if(a_size<k_size){
    //     mexErrMsgTxt("Kernel is longer than signal");
    // }
    if(n != mxGetN(prhs[1])){
        mexErrMsgTxt("Kernel and signal must have the same number of columns");
    }

    plhs[0] = mxCreateDoubleMatrix(a_size + k_size - 1, n, mxREAL);
    conv = mxGetPr(plhs[0]);

    struct threadData data[NUMTHREADS];
    pthread_t thread[NUMTHREADS];
    int tasksPerThread = (n + NUMTHREADS - 1) / NUMTHREADS;

    double *padded[NUMTHREADS];
    for(i=0; i<NUMTHREADS; i++){
        padded[i] = mxCalloc(a_size + 2 * (k_size - 1), sizeof(double));
        // padded[i] = mxCreateDoubleMatrix(1, (a_size + 2 * (k_size - 1)), mxREAL);
        data[i].padded_pt = padded[i];
        data[i].kernel = kernel;
        data[i].A = A;
        data[i].conv = conv;
        data[i].k_size = k_size;
        data[i].a_size = a_size;
        data[i].start = i*tasksPerThread;
        data[i].stop = (i+1)*tasksPerThread;
        data[i].sub = 0;
    }
    data[NUMTHREADS - 1].stop = n;

    for(i=0; i<NUMTHREADS; i++){
        pthread_create(&thread[i], NULL, convolve, &data[i]);
    }

    for(i=0; i<NUMTHREADS; i++){
        pthread_join(thread[i], NULL);
    }


    for(i=0; i<NUMTHREADS; i++){
      mxFree(padded[i]);
    }

    return;
}

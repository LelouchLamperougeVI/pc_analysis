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

#define NUMTHREADS 7

struct threadData{
    int start, stop, a_size, k_size;
    double *padded_pt, *kernel, *A, *conv;
};

void cpyColumn(double *receiver, double *giver, int col, int start, int n){ //copy matrix column to array
    for(int i=0; i<n; i++){
        memcpy(receiver + start + i, giver + i + n*col, sizeof(double));
    }
}

void *convolve(struct threadData *data){
    int start = data->start;
    int stop = data->stop;
    int a_size = data->a_size;
    int k_size = data->k_size;
    double *padded_pt = data->padded_pt;
    double *kernel = data->kernel;
    double *A = data->A;
    double *conv = data->conv;

    for(int i=start; i<stop; i++){
        cpyColumn(padded_pt, A, i, k_size - 1, a_size);
        for(int j=0; j<=(a_size + 2*(k_size-1) - k_size); j++){
            for(int k=0; k<k_size; k++){
                conv[i*(a_size+k_size-1) + j] += padded_pt[j+k] * kernel[i*k_size + k_size - k - 1];
            }
        }
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

    double *A, *kernel, *conv;
    int n, a_size, k_size;

    A = mxGetPr(prhs[0]);
    kernel = mxGetPr(prhs[1]);
    n = mxGetN(prhs[0]);
    a_size = mxGetM(prhs[0]);
    k_size = mxGetM(prhs[1]);

    if(a_size<k_size){
        mexErrMsgTxt("Kernel is longer than signal");
    }
    if(n != mxGetN(prhs[1])){
        mexErrMsgTxt("Kernel and signal must have the same number of columns");
    }

    plhs[0] = mxCreateDoubleMatrix(a_size + k_size - 1, n, mxREAL);
    conv = mxGetPr(plhs[0]);

    struct threadData data[NUMTHREADS];
    pthread_t thread[NUMTHREADS];
    int tasksPerThread = (n + NUMTHREADS - 1) / NUMTHREADS;

    mxArray *padded;
    for(int i=0; i<NUMTHREADS; i++){
        padded = mxCreateDoubleMatrix(1, (a_size + 2 * (k_size - 1)), mxREAL);
        data[i].padded_pt = mxGetPr(padded);
        data[i].kernel = kernel;
        data[i].A = A;
        data[i].conv = conv;
        data[i].k_size = k_size;
        data[i].a_size = a_size;
        data[i].start = i*tasksPerThread;
        data[i].stop = (i+1)*tasksPerThread;
    }
    data[NUMTHREADS - 1].stop = n;

    for(int i=0; i<NUMTHREADS; i++){
        pthread_create(&thread[i], NULL, convolve, &data[i]);
    }
    
    for(int i=0; i<NUMTHREADS; i++){
        pthread_join(thread[i], NULL);
    }
    
    return;
}
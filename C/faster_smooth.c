/*************************************************************************
 * MATLAB MEX ROUTINE faster_smooth.c
 *
 * C implementation of fast_smooth.m
 *
 * Direct complaints to HaoRan Chang
 ************************************************************************/

#include "string.h"
#include "mex.h"
#include "math.h"
#include "convolve.h"
#include "pthread.h"

#define NUMTHREADS 8

double *gausswin(int N, double sd){
  double *kernel = mxCalloc(N, sizeof(double));
  int i;
  double sum = 0;
  for(i=0; i<N; i++){
    kernel[i] = -((double) N - 1)/2 + i;
    kernel[i] = exp( -pow(kernel[i], 2) / 2 / pow(sd, 2));
    sum += kernel[i];
  }
  for(i=0; i<N; i++){
    kernel[i] /= sum;
  }
  return kernel;
}

double *repmat(int n, double *A, int a_size){
  int i, j;
  double *B = mxCalloc(n*a_size, sizeof(double));
  for(i=0; i<n; i++){
    for(j=0; j<a_size; j++){
      B[i*a_size+j] = A[j];
    }
  }
  mxFree(A);
  return B;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *signal, sd, *kernel, *conv, *smoothed;
    int m, n, i, j;
    mxArray *padded;

    m=mxGetM(prhs[0]);
    n=mxGetN(prhs[0]);
    signal = mxGetPr(prhs[0]);
    sd = *(mxGetPr(prhs[1]));
    int k_size = (int) ceil(10 * sd);

    kernel = gausswin(k_size, sd);
    kernel = repmat(n, kernel, k_size);

    // conv = mxGetPr(mxCreateDoubleMatrix(m + k_size - 1, n, mxREAL));

    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    smoothed = mxGetPr(plhs[0]);

    struct threadData data[NUMTHREADS];
    pthread_t thread[NUMTHREADS];
    int tasksPerThread = (n + NUMTHREADS - 1) / NUMTHREADS;
    for(i=0; i<NUMTHREADS; i++){
      padded = mxCreateDoubleMatrix(1, (m + 2 * (k_size - 1)), mxREAL);
      data[i].padded_pt = mxGetPr(padded);
      data[i].kernel = kernel;
      data[i].A = signal;
      data[i].conv = smoothed;
      data[i].k_size = k_size;
      data[i].a_size = m;
      data[i].start = i*tasksPerThread;
      data[i].stop = (i+1)*tasksPerThread;
      data[i].sub = 1;
    }
    data[NUMTHREADS - 1].stop = n;

    for(i=0; i<NUMTHREADS; i++){
        pthread_create(&thread[i], NULL, convolve, &data[i]);
    }

    for(i=0; i<NUMTHREADS; i++){
        pthread_join(thread[i], NULL);
    }

    // for(i=0; i<n; i++){
    //   for(j=0; j<m; j++){
    //     smoothed[i*m+j] = conv[i*(m+k_size-1) + j + (int) floor(k_size/2)];
    //   }
    // }

    mxFree(kernel);
    return;
}

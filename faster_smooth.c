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

double *gausswin(int N, double sd){
  double *kernel = mxCalloc(N, sizeof(double));
  int i;
  for(i=0; i<N; i++){
    kernel[i] = -((double) N - 1)/2 + i;
    kernel[i] = exp( -pow(kernel[i], 2) / 2 / pow(sd, 2));
  }
  return kernel;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *signal, sd, *kernel, *smoothed;
    int m, n;
    mxArray *padded;

    m=mxGetM(prhs[0]);
    n=mxGetN(prhs[0]);
    signal = mxGetPr(prhs[0]);
    sd = *(mxGetPr(prhs[1]));

    kernel = gausswin((int) ceil(10 * sd), sd);

    plhs[0] = mxCreateDoubleMatrix(m + sizeof(kernel) + 1, n, mxREAL);
    smoothed = mxGetPr(plhs[0]);

    struct threadData data;
    padded = mxCreateDoubleMatrix(1, (m + 2 * (sizeof(kernel) - 1)), mxREAL);
    data.padded_pt = mxGetPr(padded);
    data.kernel = kernel;
    data.A = signal;
    data.conv = smoothed;
    data.k_size = sizeof(kernel);
    data.a_size = m;
    data.start = 0;
    data.stop = n;

    convolve(&data);

    mxFree(kernel);
    return;
}

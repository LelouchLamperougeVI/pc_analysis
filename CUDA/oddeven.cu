/**
   Odd-even sort algortihm implemented on CUDA for compilation using mexcuda on MATLAB
   Returns matrix with sorted rows and the corresponding ranks
   Usage:
   [sorted, ranks] = oddeven(A);

   Tested with MATLAB 2018a, and compiled using CUDA 9.0 and MSVC++ 2015

   Written by HaoRan Chang,
   Polaris Brain Dynamics Research Group,
   Canadian Centre for Behavioural Neuroscience,
   University of Lethbridge, AB, Canada
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "matrix.h"

#define bool int
#define true 1
#define false 0

#define BLOCK_SIZE 1024

__global__
void tMatch(double const * const x, bool * const sorted, double * const results, int const * n, int const m){
        int const tot_i = blockIdx.x * blockDim.x + threadIdx.x;
        int const block_i = blockDim.x;
        int const thread_i = threadIdx.x;
        int i;
        double buff;
        bool all_sort = true;

        __shared__ double temp[100];

        for(i = 0; i < n[thread_i]; i++)
                temp[thread_i * n[thread_i] + i] = x[blockIdx.x * m + thread_i * (n[thread_i] - n[thread_i]%2) + i];

        __syncthreads();

        while(all_sort) {
                all_sort = false;
                sorted[tot_i] = false;
                for(i = 0; i < n[thread_i] - (n[thread_i] % 2); i+=2) {
                        if(temp[thread_i * n[thread_i] + i] > temp[thread_i * n[thread_i] + i + 2]) {
                                buff = temp[thread_i * n[thread_i] + i];
                                temp[thread_i * n[thread_i] + i] = temp[thread_i * n[thread_i] + i + 2];
                                temp[thread_i * n[thread_i] + i + 2] = buff;

                                sorted[tot_i] = false;
                        }
                }
                for(i = 1; i < n[thread_i]; i+=2) {
                        if(temp[thread_i * n[thread_i] + i] > temp[thread_i * n[thread_i] + i + 2]) {
                                buff = temp[thread_i * n[thread_i] + i];
                                temp[thread_i * n[thread_i] + i] = temp[thread_i * n[thread_i] + i + 2];
                                temp[thread_i * n[thread_i] + i + 2] = buff;

                                sorted[tot_i] = false;
                        }
                }
                __syncthreads();

                for(i = 0; i < block_i; i++)
                        all_sort = all_sort + sorted[blockIdx.x * blockDim.x + i];
        }

        for(i = 0; i < n[thread_i]; i++)
                results[blockIdx.x * m + thread_i * (n[thread_i] - n[thread_i]%2) + i] = temp[thread_i * n[thread_i] + i];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        if(nrhs < 1)
                mexErrMsgTxt("There needs to be at least one input you dummy -_-\"");

        double const *x;
        double *results;
        bool *sorted;
        mxGPUArray const *x_pr, *el_gpu_pr;
        mxGPUArray *sort, *res;
        mxArray *el_pr;
        mwSize const *dims;
        mwSize dimensions[2];
        int i;
        int const *el_per_thread_gpu;

        mxInitGPU();

        x_pr = mxGPUCreateFromMxArray(prhs[0]);
        x = (double const *) (mxGPUGetDataReadOnly(x_pr));

        dims = mxGPUGetDimensions(x_pr);
        int const m = dims[0];
        int const n = dims[1];

        el_pr = mxCreateNumericMatrix(BLOCK_SIZE, 1, mxINT32_CLASS, mxREAL);
        int *el_per_thread = (int *) mxGetData(el_pr);
        for(i = 0; i < m - 1; i += 2)
                el_per_thread[i/2 % BLOCK_SIZE] += 2;
        if(m/BLOCK_SIZE > 0) {
                el_per_thread[BLOCK_SIZE - 1] += m%2;
        }else{
                el_per_thread[m/2 - 1] += m%2;
        }
        el_gpu_pr = mxGPUCreateFromMxArray(el_pr);
        el_per_thread_gpu = (int const *) mxGPUGetDataReadOnly(el_gpu_pr);

        for(i=0; i<BLOCK_SIZE; i++)
                mexPrintf("%d ", el_per_thread[i]);

        res = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(x_pr), dims,
                                  mxGPUGetClassID(x_pr), mxGPUGetComplexity(x_pr),
                                  MX_GPU_DO_NOT_INITIALIZE);

        dimensions[0] = BLOCK_SIZE;
        dimensions[1] = n;
        sort = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(x_pr), dimensions,
                                   mxINT32_CLASS, mxREAL,
                                   MX_GPU_DO_NOT_INITIALIZE);

        results = (double *) (mxGPUGetData(res));
        sorted = (bool *) (mxGPUGetData(sort));

        tMatch<<<n, BLOCK_SIZE>>>(x, sorted, results, el_per_thread_gpu, m);

        plhs[0] = mxGPUCreateMxArrayOnGPU(res);

        mxGPUDestroyGPUArray(x_pr);
        mxGPUDestroyGPUArray(res);
        mxGPUDestroyGPUArray(sort);
        mxGPUDestroyGPUArray(el_gpu_pr);
        mxDestroyArray(el_pr);

        return;
}

#include "mex.h"
#include "gpu/mxGPUArray.h"

#define bool int
#define true 1
#define false 0

__global__
void tMatch(double const * const x, bool * const sorted, double * const results, int n){
  int const tot_i = blockIdx.x * blockDim.x + threadIdx.x;
  int const block_i = blockDim.x;
  int const thread_i = threadIdx.x;
  int i;
  double buff;
  bool all_sort = true;

 __shared__ double temp[100];

 for(i = 0; i < n; i++){
   temp[thread_i * n + i] = x[thread_i * n + i];
 }

 __syncthreads();

 while(all_sort){
   all_sort = false;
   sorted[thread_i] = false;
    for(i = 0; i < n; i+=2){
      if(temp[thread_i * n + i] > temp[thread_i * n + i + 2]){
        buff = temp[thread_i * n + i];
        temp[thread_i * n + i] = temp[thread_i * n + i + 2];
        temp[thread_i * n + i + 2] = buff;

        sorted[thread_i] = false;
      }
    for(i = 1; i < n; i+=2){
      if(temp[thread_i * n + i] > temp[thread_i * n + i + 2]){
        buff = temp[thread_i * n + i];
        temp[thread_i * n + i] = temp[thread_i * n + i + 2];
        temp[thread_i * n + i + 2] = buff;

        sorted[thread_i] = false;
      }
     __syncthreads();

     for(i = 0; i < block_i; i++){
       all_sort = all_sort + sorted[i];
     }
    }
  }

  for(i = 0; i < n; i++){
    results[thread_i * n + i] = temp[thread_i * n + i];
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  double const *x;
  double *sorted, *results;
  bool *sorted
  mxGPUArray const *x_pr;
  mxGPUArray *sort, *res;
  mwSize const *dims;

  mxInitGPU();

  x_pr = mxGPUCreateFromMxArray(prhs[0]);
  x = (double const *) (mxGPUGetDataReadOnly(x_pr));

  dims = mxGPUGetDimensions(x_pr);
  int const m = dims[0];
  int const n = dims[1];

  res = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(x_pr), mxGPUGetDimensions(x_pr),
                            mxGPUGetClassID(x_pr), mxGPUGetComplexity(x_pr),
                            MX_GPU_DO_NOT_INITIALIZE);

  sort = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(x_pr), mxGPUGetDimensions(x_pr),
                            mxGPUGetClassID(x_pr), mxGPUGetComplexity(x_pr),
                            MX_GPU_DO_NOT_INITIALIZE);

  results = (double *) (mxGPUGetData(res));
  sorted = (double *) (mxGPUGetData(sort));

  tMatch<<<n, m>>>(templates, match, results);

  plhs[0] = mxGPUCreateMxArrayOnGPU(res);

  mxGPUDestroyGPUArray(temp);
  mxGPUDestroyGPUArray(mat);
  mxGPUDestroyGPUArray(res);
}

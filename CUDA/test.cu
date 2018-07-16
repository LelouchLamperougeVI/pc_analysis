#include "mex.h"
#include "gpu/mxGPUArray.h"

#define bool int
#define true 1
#define false 0

__global__
void tMatch(double const * const templates, double const * const match, double * const results){
  int const tot_i = blockIdx.x * blockDim.x + threadIdx.x;
  // int const block_i = blockDim.x;
  int const thread_i = threadIdx.x;
  int i;

 __shared__ double temp[100];

 temp[thread_i] = match[tot_i];

 __syncthreads();

  for(i = 0; i < 10; i++){
    temp[thread_i] += templates[thread_i];
   __syncthreads();
  }

  results[tot_i] = temp[thread_i];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  double const *templates;
  double const *match;
  double *results;
  mxGPUArray const *temp;
  mxGPUArray const *mat;
  mxGPUArray *res;
  mwSize const *dims;

  mxInitGPU();

  temp = mxGPUCreateFromMxArray(prhs[0]);
  templates = (double const *) (mxGPUGetDataReadOnly(temp));
  mat = mxGPUCreateFromMxArray(prhs[1]);
  match = (double const *) (mxGPUGetDataReadOnly(mat));

  dims = mxGPUGetDimensions(mat);
  int const m = dims[0];
  int const n = dims[1];

  res = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(mat), mxGPUGetDimensions(mat), mxGPUGetClassID(mat), mxGPUGetComplexity(mat), MX_GPU_DO_NOT_INITIALIZE);
  results = (double *) (mxGPUGetData(res));

  tMatch<<<n, m>>>(templates, match, results);

  plhs[0] = mxGPUCreateMxArrayOnGPU(res);

  mxGPUDestroyGPUArray(temp);
  mxGPUDestroyGPUArray(mat);
  mxGPUDestroyGPUArray(res);
}

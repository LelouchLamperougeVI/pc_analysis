/* The default MATLAB function xcorr2 does not allow you specify a maxlag range
 * This could result in very poor performance for certain applications that do not require convolving the entire matrices
 * This algorithm does not make use of the Convolution Theorem. Therefore, for large convolutions, use built-in function instead.
 *
 * Originally written for an old Intel Westmere architecture (dual socket Xeon E5620) using SSE4.1
 * Performance can be improved by reimplementing for AVX-512
 *
 * Compilation:
 *    mex -R2018a CXXFLAGS='$CXXGLAGS -msse4.1' matxcorr2.cpp corr.cpp
 *
 * Usage:
 *    r = matxcorr2(A, B, int16(maxlag), numthreads)
 *
 * Inputs:
 *    A, B            input matrices; B is moving
 *    maxlag          a 1x2 vector specifying the lags for each dimension
 *    numthreads      number of compute threads - note: this function does leverage hyperthreading :)
 *
 * Outputs:
 *    r               raw correlation coefficients without any sort of normalization
 *
 * To be used in conjunction with the accompanying matlab function bxcorr2().
 *
 * By HaoRan Chang, Ph.D. candidate
 * Canadian Centre for Behavioural Neuroscience,
 * Department of Neuroscience, University of Lethbridge,
 * Lethbridge, Alberta, Canada
 *
 * 2019
 *
 */
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include <matrix.h>
#include "corr.h"

namespace corr {
  double *A, *B, *ret;
  int16_t *maxlag;
  int dims[4], S[2], *shift, mod[2];
  std::mutex mute;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if(nrhs < 4)
    mexErrMsgTxt("matxcorr2 requires 4 inputs");
  if(nrhs > 4)
    mexWarnMsgTxt("ignoring extra input arguments");
  if(nlhs > 3)
    mexErrMsgTxt("too many output arguments");

  if( !mxIsDouble(prhs[0]) | !mxIsDouble(prhs[1]) | !mxIsDouble(prhs[3]) )
    mexErrMsgTxt("input matrices A and B, and numThreads must be of type double");
  if( !mxIsInt16(prhs[2]) )
    mexErrMsgTxt("maxlag must be given as int16_t");

  int i, numThreads;

  using namespace corr;
  using namespace std;

  A = mxGetDoubles(prhs[0]);
  B = mxGetDoubles(prhs[1]);
  if(mxGetNumberOfDimensions(prhs[0]) != 2 | mxGetNumberOfDimensions(prhs[1]) != 2)
    mexErrMsgTxt("the input matrices must contain 2 and only 2 dimensions");
  dims[0] = (int) mxGetM(prhs[0]);
  dims[1] = (int) mxGetN(prhs[0]);
  dims[2] = (int) mxGetM(prhs[1]);
  dims[3] = (int) mxGetN(prhs[1]);

  S[0] = dims[0] + dims[2] - 1;
  S[1] = dims[1] + dims[3] - 1;

  maxlag = mxGetInt16s(prhs[2]);
  numThreads = (mwSize) mxGetScalar(prhs[3]);

  if((S[0]/2 - (int) maxlag[0]) < 0 | (S[1]/2 - (int) maxlag[1]) < 0)
    mexErrMsgTxt("maxlag exceeds M+N-1 maximum allowence");

  mod[0] = 2 * (int) maxlag[0] + S[0] % 2;
  mod[1] = 2 * (int) maxlag[1] + S[1] % 2;
  if(!mod[0] || !mod[1])
    mexErrMsgTxt("behaviour is undefined for maxlag = 0 with even-odd matrix-pairs");

  plhs[0] = mxCreateDoubleMatrix(mod[0], mod[1], mxREAL);
  ret = mxGetDoubles(plhs[0]);

  shift = (int *) calloc(2, sizeof(int));
  thread workers[numThreads];
  for(i=0; i<numThreads; i++)
    workers[i] = thread(fWorker);
  for(i=0; i<numThreads; i++)
    workers[i].join();

  }

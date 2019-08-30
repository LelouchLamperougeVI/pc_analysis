/* The default MATLAB function `circshift' doesn't allow you to define a vector to shift individual rows/columns.
 * Not that big of a deal but still a pain in the ass and can be very slow for permutation test.
 *
 * Compilation:
 *    mex -R2018a matcircshift.cpp
 *
 * Usage:
 *    Y = matcircshift(A, k, dim)
 *
 * Inputs:
 *    A:    matrix to shift
 *    k:    vector of shift amount; has to be the same length as the colum dimension
 *
 * Output:
 *    Y: shifted matrix
 *
 * To be used in conjunction with the accompanying matlab function bcircshift().
 *
 *
 * By HaoRan Chang, Ph.D. candidate
 * Canadian Centre for Behavioural Neuroscience,
 * Department of Neuroscience, University of Lethbridge,
 * Lethbridge, Alberta, Canada
 *
 * 2019
 *
 */

#include "mex.h"
#include <thread>
#include <mutex>

namespace shifter {
  unsigned int numThreads = std::thread::hardware_concurrency(); //automatically detect number of threads available, only tested on linux
  int M, N;
  double *A, *k, *Y;
  std::mutex mute;

  class shifter;
  void shift();
};

class shifter::shifter {
  public:
    void operator()() {
      do {
        mute.lock();
        int index = current;
        current++;
        mute.unlock();
        if(index >= N)
          return;

        for(int i = 0; i < M; i++)
          Y[i + M*index] = A[(i + (int) k[index]) % M + M*index];
      } while(1);
    }
  private:
    int current;
};

void shifter::shift() {
  shifter queue;
  std::thread workers[numThreads];
  for(int i = 0; i < numThreads; i++)
    workers[i] = std::thread(queue);
  for(int i = 0; i < numThreads; i++)
    workers[i].join();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if(nrhs != 2)
    mexErrMsgTxt("two function arguments must be passed");
  if(nlhs > 1)
    mexErrMsgTxt("too many output arguments");
  if( !mxIsDouble(prhs[0]) | !mxIsDouble(prhs[1]) )
    mexErrMsgTxt("input arguments must be of type `double'");

  shifter::A = mxGetDoubles(prhs[0]);
  shifter::k = mxGetDoubles(prhs[1]);
  shifter::M = mxGetM(prhs[0]);
  shifter::N = mxGetN(prhs[0]);

  if((int) mxGetNumberOfElements(prhs[1]) == 1) {
    double temp = *shifter::k;
    shifter::k = new double[shifter::N];
    for(int i = 0; i < shifter::N; i++) shifter::k[i] = temp;
  }else{
    if((int) mxGetNumberOfElements(prhs[1]) != shifter::N)
      mexErrMsgTxt("input vector `k' does not match the length of the columns dimension");
  }

  plhs[0] = mxCreateDoubleMatrix(shifter::M, shifter::N, mxREAL);
  shifter::Y = mxGetDoubles(plhs[0]);

  shifter::shift();

  if((int) mxGetNumberOfElements(prhs[1]) == 1)
    delete [] shifter::k;
}

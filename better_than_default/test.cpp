#include "corr.h"
#include <stdlib.h>

namespace corr {
  double *A, *B, *ret;
  int16_t *maxlag;
  int dims[4], S[2], *shift, mod[2];
  std::mutex mute;
}

int main() {
  int i, M=100, N=100;

  using namespace corr;
  using namespace std;
  A = (double *) aligned_alloc(16, M * N * sizeof(double));
  B = (double *) aligned_alloc(16, M * N * sizeof(double));

  for(i=0; i<M*N; i++){
    A[i]=1;
    B[i]=1;
  }

  dims[0] = M;
  dims[1] = N;
  dims[2] = M;
  dims[3] = N;

  S[0] = dims[0] + dims[2] - 1;
  S[1] = dims[1] + dims[3] - 1;

  maxlag = (int16_t *) malloc(2 * sizeof(int16_t));
  maxlag[0] = 99;
  maxlag[1] = 99;

  int numThreads = 16;

  ret = (double *) calloc(S[0] * S[1], sizeof(double));

  mod[0] = 2 * (int) maxlag[0] + S[0] % 2;
  mod[1] = 2 * (int) maxlag[1] + S[1] % 2;

  shift = (int *) calloc(2, sizeof(int));
  thread workers[numThreads];
  for(i=0; i<numThreads; i++)
    workers[i] = thread(fWorker);
  for(i=0; i<numThreads; i++)
    workers[i].join();
}

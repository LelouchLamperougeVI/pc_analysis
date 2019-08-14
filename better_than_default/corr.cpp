#include "corr.h"

void correlate(int *index) {
  /***
  algo:
    A: arg_max(shift - B + 1, 0) --> arg_min(shift, A - 1)
    B: arg_max(S - shift - A, 0) --> arg_min(S - shift - 1, B - 1)
  ***/

  using namespace std;
  using namespace corr;

  int i, j, offset[4], cap[2];

  double *ret_pt = ret + index[0] + index[1]*mod[0];
  index[0] += S[0]/2 - (int) maxlag[0];
  index[1] += S[1]/2 - (int) maxlag[1];

  offset[0] = max(index[0] - dims[2] + 1, 0);
  offset[1] = max(S[0] - index[0] - dims[0], 0);
  offset[2] = max(index[1] - dims[3] + 1, 0);
  offset[3] = max(S[1] - index[1] - dims[1], 0);

  cap[0] = min(index[0], dims[0] - 1) - offset[0] + 1;
  cap[1] = min(index[1], dims[1] - 1) - offset[2] + 1;

  double *A_pt = A + offset[0] + offset[2]*dims[0];
  double *B_pt = A + offset[1] + offset[3]*dims[2];

  __m128d A_store;
  __m128d B_store;

  __m128d buff = _mm_setzero_pd();
  for(j=0; j<cap[1]; j++){
    for(i=0; i<cap[0]; i+=XMM_SIZE){
      A_store = _mm_loadu_pd(A_pt+i);
      B_store = _mm_loadu_pd(B_pt+i);
      buff += _mm_dp_pd( A_store, B_store, 0xF1);
    }
    if(cap[0] % XMM_SIZE)
      *ret += *(A_pt+i-1) * *(B_pt+i-1);
    A_pt = A + offset[0] + (j+offset[2])*dims[0];
    B_pt = A + offset[1] + (j+offset[3])*dims[2];
  }
  *ret_pt += _mm_cvtsd_f64(buff);
}


void fWorker() {
  using namespace corr;

  int offset[2];
  std::function<void ()> shifter = [](){
    shift[1] += ++shift[0] % mod[0] ? 0 : !( shift[0] = 0 );
  };
  do{
    mute.lock();
    if( shift[1] == mod[1] ) {
      mute.unlock();
      return;
    }
    offset[0] = shift[0];
    offset[1] = shift[1];
    shifter();
    mute.unlock();
    correlate(offset);
  } while(1);
}


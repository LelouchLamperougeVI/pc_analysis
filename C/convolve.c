#include "string.h"
#include "convolve.h"

void cpyColumn(double *receiver, double *giver, int col, int start, int n){
    int i;
    for(i=0; i<n; i++){
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
    int sub = data->sub;

    int i, j, k;

    if(sub){
      for(i=start; i<stop; i++){
          cpyColumn(padded_pt, A, i, k_size - 1 - k_size/2, a_size);
          for(j=0; j<a_size; j++){
              for(k=0; k<k_size; k++){
                  conv[i*(a_size) + j] += padded_pt[j+k] * kernel[i*k_size + k_size - k - 1];
              }
          }
      }
    }else{
      for(i=start; i<stop; i++){
          cpyColumn(padded_pt, A, i, k_size - 1, a_size);
          for(j=0; j<(a_size + 2*(k_size-1) - k_size + 1); j++){
              for(k=0; k<k_size; k++){
                  conv[i*(a_size+k_size-1) + j] += padded_pt[j+k] * kernel[i*k_size + k_size - k - 1];
              }
          }
      }
    }
}

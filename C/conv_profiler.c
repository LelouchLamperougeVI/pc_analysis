#include <stdio.h>
#include <math.h>
#include <string.h>
#include "convolve.h"

void *main(){
  struct threadData data;
  double *signal, *kernel, *conv, *pad;
  int a_size, k_size, i, j;

  int n=100;
  a_size=10000;
  k_size=1000;

  signal = calloc(a_size*n, sizeof(double));
  kernel = calloc(k_size*n, sizeof(double));
  pad = calloc(a_size+(k_size-1)*2, sizeof(double));
  conv = calloc(n*(a_size+k_size-1), sizeof(double));

  for(i=0; i<a_size*n; i++){
    signal[i] = (double) i;
  }
  for(i=0; i<k_size*n; i++){
    kernel[i] = (double) i;
  }

  data.padded_pt = pad;
  data.kernel = kernel;
  data.A = signal;
  data.conv = conv;
  data.k_size = k_size;
  data.a_size = a_size;
  data.start = 0;
  data.stop = n;

  convolve(&data);

  for(i=0; i<n; i++){
    for(j=0; j<(a_size+k_size-1); j++){
      printf("%d \t", (int) conv[i*(a_size+k_size-1) + j]);
    }
    printf("\n");
  }

  free(signal); free(kernel); free(pad); free(conv);

  return;
}

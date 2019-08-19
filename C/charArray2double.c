/*************************************************************************
 * MATLAB MEX ROUTINE charArray2double.c
 *
 * Direct complaints to HaoRan Chang
 ************************************************************************/

#include "string.h"
#include "stdlib.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int m, n, i;
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
  double *d_array = mxGetPr(plhs[0]);

  char *C = mxMalloc(m*n+1);
  mxGetString(prhs[0], C, m*n+1);

  char *buff = mxMalloc(m+1);

  for(i = 0; i < n; i++){
    memcpy(buff, C + i*m, m+1);
    d_array[i] = strtod(buff, NULL);
  }

  mxFree(C); mxFree(buff);

  return;
}

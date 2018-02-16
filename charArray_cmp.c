/*************************************************************************
 * MATLAB MEX ROUTINE charArray_cmp.c
 *
 * Direct complaints to HaoRan Chang
 ************************************************************************/

#include "string.h"
#include "stdlib.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int m, n, s, i;
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  s = mxGetNumberOfElements(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
  double *d_array = mxGetPr(plhs[0]);

  char *C = mxMalloc(m*n*sizeof(char)+1);
  mxGetString(prhs[0], C, m*n*sizeof(char)+1);
  char *target = mxMalloc(s*sizeof(char)+1);
  mxGetString(prhs[1], target, s*sizeof(char)+1);

  char *buff = mxCalloc(m, sizeof(char));
  char *trail;

  for(i = 0; i < n; i++){
    memcpy(buff, C + i*m, m*sizeof(char));
    trail = buff + strlen(buff) - 1;
    while(isspace((unsigned char) *trail)) trail--;
    *(trail+1) = 0;
    if(strcmp(buff, target)){
      d_array[i] = 0;
    } else{
      d_array[i] = 1;
    }
  }

  mxFree(C); mxFree(buff); mxFree(target);

  return;
}

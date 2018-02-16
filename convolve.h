#ifndef _CONV
#define _CONV

#include "string.h"
#include "mex.h"

struct threadData{
    int start, stop, a_size, k_size;
    double *padded_pt, *kernel, *A, *conv;
};

void *convolve(struct threadData *data);

#endif

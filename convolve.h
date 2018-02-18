#ifndef _CONV
#define _CONV

#include "string.h"

struct threadData{
    int start, stop, a_size, k_size;
    double *padded_pt, *kernel, *A, *conv;
    int sub; // only central part of convolution
};

void *convolve(struct threadData *data);

#endif

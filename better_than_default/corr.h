/*
 * The current implementation requires SSE4.1
 * Honestly, I don't even know why I would mention this; I'm probably the only person who's still using such old ISA.
 * Note to self: If you ever upgrade to a newer processor, optimize the code with AVX-512.
 */

#ifndef CORR_H
#define CORR_H

#define XMM_SIZE 2

#include <algorithm>
#include <functional>
#include <thread>
#include <mutex>
#include <smmintrin.h>

namespace corr {
  extern double *A, *B, *ret;
  extern int16_t *maxlag;
  extern int dims[4], S[2], *shift, mod[2];

  extern std::mutex mute;
}

void correlate(int *index);

void fWorker();

#endif

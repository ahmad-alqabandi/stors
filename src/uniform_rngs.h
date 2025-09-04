
#include "R_stors.h"

#ifndef UNIFORM_RNGS_H
#define UNIFORM_RNGS_H

// KISS RNG macros

// Global RNG state
static struct {
  uint32_t jcong;
  uint32_t jsr;
  uint32_t w;
  uint32_t z;
} rng = {380116160, 123456789, 521288629, 362436069};

//Define the KISS RNG macros
#define UINT32_MAX_DOUBLE 1.0 / (double)UINT32_MAX

#define znew (rng.z = 36969 * (rng.z & 65535) + (rng.z >> 16))
#define wnew (rng.w = 18000 * (rng.w & 65535) + (rng.w >> 16))
#define MWC  ((znew << 16) + wnew)
#define SHR3 (jz = rng.jsr, rng.jsr ^= (rng.jsr << 13), rng.jsr ^= (rng.jsr >> 17), rng.jsr ^= (rng.jsr << 5), jz + rng.jsr)
#define CONG (rng.jcong = 69069 * rng.jcong + 1234567)
#define KISS ((MWC ^ CONG) + SHR3)
#define KISS_UNIF (KISS * UINT32_MAX_DOUBLE)

// Function declarations
void KISS_set_seed(uint32_t seed);
void KISS_set_seed_R(int *seed);
void KISS_unif(double *results, int *n);



// INIT Uniform RNG function pointer
double (*u_rng)();


#endif

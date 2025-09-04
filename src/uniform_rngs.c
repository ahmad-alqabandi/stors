#include "uniform_rngs.h"


// KISS RNG functions

uint32_t jz;

void KISS_set_seed(uint32_t seed) {
  rng.jsr = 123456789;
  if (rng.jsr != seed) {
    rng.jsr ^= seed;
  }
  rng.z     = 362436069;
  rng.w     = 521288629;
  rng.jcong = 380116160;
}

void KISS_set_seed_R(int *seed) {
  KISS_set_seed((uint32_t)*seed);
}

void KISS_unif(double *results, int *n) {
  for (int i = 0; i < *n; i++) {
    results[i] = KISS_UNIF;
  }
}

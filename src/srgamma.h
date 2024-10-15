#include "R_stors.h"


#ifndef SRGAMMA_H
#define SRGAMMA_H


SEXP srgamma(SEXP s_size);
SEXP srgamma_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP srgamma_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir);
  
#endif



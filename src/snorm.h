#include "R_stors.h"


#ifndef SNORM_H
#define SNORM_H


SEXP srnorm(SEXP s_size, SEXP Rpassed_params);
SEXP srnorm_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP srnorm_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir);
  
#endif



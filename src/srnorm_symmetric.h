#include "R_stors.h"


#ifndef SNORM_S_H
#define SNORM_S_H

SEXP srnorm_symmetric(SEXP s_size);
SEXP srnorm_symmetric_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP srnorm_symmetric_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir);
  
#endif



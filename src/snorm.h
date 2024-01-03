#include "R_stors.h"


#ifndef SNORM_H
#define SNORM_H


SEXP srnorm(SEXP s_size);
SEXP srnorm_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP srnorm_trunc(SEXP s_size, SEXP l, SEXP u);
  
#endif



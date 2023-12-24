#ifndef SNORM_H
#define SNORM_H


#include <R.h>
#include <Rinternals.h>

#include "cache.h"

#define R_RETURN_NULL return(R_NilValue);


SEXP srnorm(SEXP s_size);
SEXP srnorm_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP srnorm_trunc(SEXP s_size, SEXP l, SEXP u);
  
#endif



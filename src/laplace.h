#ifndef LAPLACE_H
#define LAPLACE_H


#include <R.h>
#include <Rinternals.h>

#include "cache.h"

#define R_RETURN_NULL return(R_NilValue);


SEXP laplace(SEXP s_size);
SEXP laplace_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP laplace_trunc(SEXP s_size, SEXP Rl, SEXP Rr);
  
#endif

#ifndef SNORM_H
#define SNORM_H


#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP srnorm(SEXP s_size);
SEXP slaplace(SEXP s_size);
SEXP rLaplace_c(SEXP s_size);
  
#endif

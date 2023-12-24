
#ifndef STORS_H
#define STORS_H

#include <R.h>
#include <Rinternals.h>

#include "cache.h"

#define R_RETURN_NULL return(R_NilValue);

SEXP print_cached_grids(void);
SEXP stors(SEXP s_size, SEXP R_Cnum, SEXP Rf, SEXP Renv);
  
#endif

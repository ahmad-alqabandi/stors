
#ifndef STORS_H
#define STORS_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP print_cached_grids(void);
SEXP stors(SEXP s_size, SEXP R_Cnum, SEXP Rf, SEXP Renv);
  
#endif


#ifndef STORS_H
#define STORS_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP print_cached_grids(void);
SEXP stors(SEXP s_size, SEXP Rx, SEXP Rs_upper_lower, SEXP Rp_a, SEXP Rm, SEXP Rnormalized_areas, SEXP Runif_s, SEXP Rs_upper, SEXP Rlts, SEXP Rrts, SEXP Rf, SEXP Renv);

#endif

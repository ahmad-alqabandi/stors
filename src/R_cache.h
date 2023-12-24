#ifndef RCACH_H
#define RCACH_H

#include <R.h>
#include <Rinternals.h>

#include "cache.h"

#define R_RETURN_NULL return(R_NilValue);

SEXP cache_grid(SEXP R_Cnum, SEXP R_x, SEXP R_s_upper, SEXP R_p_a, SEXP R_s_upper_lower, SEXP R_areas, SEXP R_steps_number, SEXP R_sampling_probabilities, SEXP R_unif_scaler, SEXP R_lt_properties, SEXP R_rt_properties, SEXP Ralpha);
SEXP free_cache(void);

#endif

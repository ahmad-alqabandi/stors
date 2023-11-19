#ifndef RCACH_H
#define RCACH_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP cache_grid(SEXP R_Cnum, SEXP R_x, SEXP R_S_upper, SEXP R_p_a, SEXP R_s_upper_lower, SEXP R_areas, SEXP R_strps_number, SEXP R_sampling_probabilities, SEXP R_unif_scaler, SEXP R_lt_properties, SEXP R_rt_properties);
SEXP free_cache(void);

#endif

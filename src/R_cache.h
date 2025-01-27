#include "R_stors.h"

#ifndef RCACH_H
#define RCACH_H

#include "cache.h"

SEXP cache_grid(SEXP R_Cnum, SEXP R_x, SEXP R_s_upper,
                SEXP R_p_a, SEXP R_s_upper_lower, SEXP R_areas,
                SEXP R_steps_number, SEXP R_sampling_probabilities,
                SEXP R_unif_scaler, SEXP R_lt_properties, SEXP R_rt_properties,
                SEXP Ralpha, SEXP Rsymmetric, SEXP Rparams, SEXP Rn_params,
                SEXP Rlb, SEXP lower);
SEXP free_cache(void);
SEXP free_cache_cnum( SEXP Rcnum);
  
#endif

#include "R_stors.h"


#ifndef SNORM_H
#define SNORM_H

SEXP srnorm_scaled_check(SEXP s_size, SEXP Rpassed_params, SEXP Rresults );

SEXP srnorm_scaled(SEXP s_size, SEXP Rpassed_params);
SEXP srnorm_sym_scaled(SEXP s_size, SEXP Rpassed_params);

SEXP srnorm_scaled_inplace( SEXP Rpassed_params, SEXP Rresults);
SEXP srnorm_sym_scaled_inplace( SEXP Rpassed_params, SEXP Rresults);




SEXP srnorm_custom_check(SEXP s_size, SEXP Rresults);

SEXP srnorm_custom(SEXP s_size);
SEXP srnorm_sym_custom(SEXP s_size);

SEXP srnorm_custom_inplace(SEXP Rresults);
SEXP srnorm_sym_custom_inplace(SEXP Rresults);


#endif



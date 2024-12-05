
#include "R_stors.h"

#ifndef LAPLACE_H
#define LAPLACE_H

SEXP srlaplace_scaled_check(SEXP s_size, SEXP Rpassed_params, SEXP Rresults );

SEXP srlaplace_scaled(SEXP s_size, SEXP Rpassed_params);
SEXP srlaplace_sym_scaled(SEXP s_size, SEXP Rpassed_params);

SEXP srlaplace_scaled_inplace( SEXP Rpassed_params, SEXP Rresults);
SEXP srlaplace_sym_scaled_inplace( SEXP Rpassed_params, SEXP Rresults);




SEXP srlaplace_custom_check(SEXP s_size, SEXP Rresults);

SEXP srlaplace_custom(SEXP s_size);
SEXP srlaplace_sym_custom(SEXP s_size);

SEXP srlaplace_custom_inplace(SEXP Rresults);
SEXP srlaplace_sym_custom_inplace(SEXP Rresults);


#endif

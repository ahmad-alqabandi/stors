#include "R_stors.h"


#ifndef SEXP_H
#define SEXP_H


SEXP srexp_scaled(SEXP s_size, SEXP Rpassed_params);
SEXP srexp_scaled_inplace( SEXP Rpassed_params, SEXP Rresults);
SEXP srexp_scaled_check(SEXP s_size, SEXP Rpassed_params, SEXP Rresults );


SEXP srexp_custom(SEXP s_size);
SEXP srexp_custom_inplace(SEXP Rresults);
SEXP srexp_custom_check(SEXP s_size, SEXP Rresults);


#endif



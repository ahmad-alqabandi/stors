
#include "R_stors.h"

#ifndef LAPLACE_H
#define LAPLACE_H

SEXP srlaplace_scaled_check(SEXP s_size, SEXP Rpassed_params );

SEXP srlaplace_scaled(SEXP s_size, SEXP Rpassed_params);
SEXP srlaplace_sym_scaled(SEXP s_size, SEXP Rpassed_params);




SEXP srlaplace_custom_check(SEXP s_size);

SEXP srlaplace_custom(SEXP s_size);
SEXP srlaplace_sym_custom(SEXP s_size);





SEXP srlaplace_trunc_nav(SEXP Rlx, SEXP Rrx, SEXP Rgrid_number);

SEXP srlaplace_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir, SEXP Rgrid_number);

#endif

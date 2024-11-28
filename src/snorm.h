#include "R_stors.h"


#ifndef SNORM_H
#define SNORM_H





SEXP srnorm_scaled_check(SEXP s_size, SEXP Rpassed_params );

SEXP srnorm_scaled(SEXP s_size, SEXP Rpassed_params);
SEXP srnorm_sym_scaled(SEXP s_size, SEXP Rpassed_params);




SEXP srnorm_custom_check(SEXP s_size);

SEXP srnorm_custom(SEXP s_size);
SEXP srnorm_sym_custom(SEXP s_size);





SEXP srnorm_trunc_nav(SEXP Rlx, SEXP Rrx, SEXP Rgrid_number);

SEXP srnorm_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir, SEXP Rgrid_number);
  
#endif



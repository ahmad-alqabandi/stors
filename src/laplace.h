
#include "R_stors.h"

#ifndef LAPLACE_H
#define LAPLACE_H

SEXP laplace(SEXP s_size, SEXP Rpassed_params);
SEXP laplace_trunc_nav(SEXP Rlx, SEXP Rrx, SEXP Rgrid_number);
SEXP laplace_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr, SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir, SEXP Rgrid_number);

#endif


#include "R_stors.h"

#ifndef LAPLACE_H
#define LAPLACE_H

SEXP laplace(SEXP s_size);
SEXP laplace_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP laplace_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir);

#endif

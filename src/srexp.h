#include "R_stors.h"


#ifndef SEXP_H
#define SEXP_H

SEXP srexp(SEXP s_size);
SEXP srexp_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP srexp_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir);

#endif



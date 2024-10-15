#include "R_stors.h"


#ifndef SRCHISQ_H
#define SRCHISQ_H


SEXP srchisq(SEXP s_size);
SEXP srchisq_trunc_nav(SEXP Rlx, SEXP Rrx);
SEXP srchisq_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir);
  
#endif



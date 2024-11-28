#include "R_stors.h"


#ifndef SRCHISQ_H
#define SRCHISQ_H


SEXP srchisq(SEXP s_size, SEXP Rpassed_params);
SEXP srchisq_trunc_nav(SEXP Rlx, SEXP Rrx, SEXP Rgrid_number);
SEXP srchisq_trunc(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir, SEXP Rgrid_number);
  
#endif



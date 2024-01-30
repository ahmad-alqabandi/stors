
#include "R_stors.h"

#ifndef STORS_H
#define STORS_H

SEXP stors(SEXP s_size, SEXP R_Cnum, SEXP Rf, SEXP Renv);
SEXP stors_trunc_nav(SEXP Rcnum ,SEXP Rlx, SEXP Rrx);
SEXP stors_trunc(SEXP s_size, SEXP Rcnum, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir,SEXP Rf, SEXP Renv);
  
#endif

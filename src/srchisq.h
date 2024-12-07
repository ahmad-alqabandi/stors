#include "R_stors.h"


#ifndef SRCHISQ_H
#define SRCHISQ_H

SEXP srchisq_custom_check(SEXP s_size, SEXP Rresults);

SEXP srchisq_custom(SEXP s_size);

SEXP srchisq_custom_inplace(SEXP Rresults);

#endif



#include "R_stors.h"


#ifndef SRPARETO_H
#define SRPARETO_H


SEXP srpareto_custom_check(SEXP s_size, SEXP Rresults);

SEXP srpareto_custom(SEXP s_size);

SEXP srpareto_custom_inplace(SEXP Rresults);

#endif

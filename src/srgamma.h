#include "R_stors.h"


#ifndef SRGAMMA_H
#define SRGAMMA_H


SEXP srgamma_custom_check(SEXP s_size, SEXP Rresults);

SEXP srgamma_custom(SEXP s_size);

SEXP srgamma_custom_inplace(SEXP Rresults);

#endif

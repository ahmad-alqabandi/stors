#include "R_stors.h"


#ifndef SRBETA_H
#define SRBETA_H

SEXP srbeta_custom_check(SEXP s_size, SEXP Rresults);

SEXP srbeta_custom(SEXP s_size);

SEXP srbeta_custom_inplace(SEXP Rresults);

#endif

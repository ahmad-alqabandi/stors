#include "macro_var.h"


#include "R_stors.h"

#ifndef TABL_H
#define TABL_H


#define CNUM 1

#define SCALE(x) x * pp[1] + pp[0]

#define L_TAIL ARS

#define R_TAIL ARS

#define F(x) ( ( g->params[1]  ) * exp(-0.5 * ((x - g->params[0]) * g->params[1]) * ((x - g->params[0]) * g->params[1])))


SEXP tabl_srnorm_custom_inplace(SEXP Rresults);


#endif

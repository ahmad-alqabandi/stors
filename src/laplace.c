#include "laplace.h"
#include "macrodef.h"

#define CNUM 1

#define L_TAIL IT

#define R_TAIL IT

#define L_ITF(u)(log(2*u))

#define R_ITF(u)(-log(2 - 2 * u))

#define F(sample)(0.5 * exp(-fabs(sample)))

SEXP laplace(SEXP s_size){
#include "stors_body.h"
}

#undef CNUM

#undef L_TAIL

#undef R_TAIL
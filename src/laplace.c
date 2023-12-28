#include "laplace.h"
#include "macrodef.h"

#define CNUM 1

#define L_TAIL IT

#define R_TAIL IT

#define L_ITF(u)(log(2*u))

#define R_ITF(u)(-log(2 - 2 * u))

#define F(x)(0.5 * exp(-fabs(x)))

#define CDF(x) ( (x <= 0) ? (0.5 * exp(x)) : (1 - 0.5 * exp(-x)) )

SEXP laplace(SEXP s_size){
#include "stors_body.h"
}

// =================================

SEXP laplace_trunc_nav(SEXP Rlx, SEXP Rrx){
#include "trunc_stors_body.h"
}

// =================================

#define TRUNC

SEXP laplace_trunc(SEXP s_size, SEXP Rl, SEXP Rr)
{
  
#include "stors_body.h"
  
}

#undef TRUNC

// =================================


#undef CNUM

#undef L_TAIL

#undef R_TAIL
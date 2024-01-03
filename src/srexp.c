#include "srexp.h"
#include "macro_var.h"

#define NAME srexp

#define CNUM 2

#define R_TAIL IT

#define R_ITF(u)(1-exp(-u))

#define F(x)(exp(-x))

// =================================

#include "stors_sample.c"

// =================================

SEXP srexp_trunc_nav(SEXP Rlx, SEXP Rrx){
#include "trunc_stors_body.h"
}

// =================================

#define TRUNC

#include "stors_sample.c"

#undef TRUNC

// =================================


#undef CNUM

#undef R_TAIL







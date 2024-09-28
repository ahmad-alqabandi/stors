#include "laplace.h"
#include "macro_var.h"

#define NAME laplace

#define CNUM 2

#define L_TAIL IT

#define R_TAIL IT

#define L_ITF(u)(0.5 * exp(u))

#define R_ITF(u)(1 - (0.5 * exp(-u)))

#define F(x)(0.5 * exp(-fabs(x)))

#define CDF(x) ( (x <= 0) ? (0.5 * exp(x)) : (1 - 0.5 * exp(-x)) )

#include "stors_sample.c"

// =================================

#include "stors_trunc_nav.c"

// =================================

#include "stors_sample_trunc.c"

// =================================


#undef CNUM

#undef L_TAIL

#undef R_TAIL
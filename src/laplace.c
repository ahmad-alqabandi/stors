#include "laplace.h"
#include "macro_var.h"

#define NAME laplace

#define CNUM 2

#define L_TAIL IT

#define R_TAIL IT

#define SYMMETRIC FALSE

#define L_ITF(u) (g.params[0] + g.params[1] * log(2 * (u)))

#define R_ITF(u) (g.params[0] - g.params[1] * log(2 - 2 * (u)))

#define F(x) (1.0 / (2.0 * g.params[1]) * exp(-fabs((x) - (g.params[0])) / (g.params[1])))

#define CDF(x) (((x) <= (g.params[0])) ? (0.5 * exp((x - g.params[0]) / g.params[1])) : (1 - 0.5 * exp(-(x - g.params[0]) / g.params[1])))

#include "stors_sample.c"

// =================================

#include "stors_trunc_nav.c"

// =================================

#include "stors_sample_trunc.c"

// =================================


#undef CNUM

#undef L_TAIL

#undef R_TAIL

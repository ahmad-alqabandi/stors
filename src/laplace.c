#include "laplace.h"
#include "macro_var.h"

#define NAME laplace

#define CNUM 3

#define CNUM_SCALED 4

#define L_TAIL IT

#define R_TAIL IT

#define L_ITF(u) (g->params[0] + g->params[1] * log(2 * (u)))

#define R_ITF(u) (g->params[0] - g->params[1] * log(2 - 2 * (u)))

#define F(x) (1.0 / (2.0 * g->params[1]) * exp(-fabs((x) - (g->params[0])) / (g->params[1])))

#define CDF(x) (((x) <= (g->params[0])) ? (0.5 * exp((x - g->params[0]) / g->params[1])) : (1 - 0.5 * exp(-(x - g->params[0]) / g->params[1])))

#define FLIP_SAMPLE(sample, flip) flip ? g->symmetric - (sample - g->symmetric) : sample

#define SCALABLE

#define SCALE(sample) sample * pp[1] + pp[0]

#include "stors_sample.c"

// =================================

#include "stors_trunc_nav.c"

// =================================

#include "stors_sample_trunc.c"

// =================================


#undef CNUM

#undef L_TAIL

#undef R_TAIL

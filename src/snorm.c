#include "snorm.h"
#include "macro_var.h"

#define NAME srnorm

#define CNUM 1

#define CNUM_SCALED 2

#define L_TAIL ARS

#define R_TAIL ARS

#define FLIP_SAMPLE(sample, flip) flip ? g.symmetric - (sample - g.symmetric) : sample

#define SCALABLE

#define SCALE(sample) sample * pp[1] + pp[0]

#define F(sample) (1.0 / (g.params[1] * 2.50662827463) * exp(-0.5 * ((sample - g.params[0]) / g.params[1]) * ((sample - g.params[0]) / g.params[1])))

// =================================

#include "stors_sample.c"

// =================================

#include "stors_trunc_nav.c"

// =================================

#include "stors_sample_trunc.c"

// =================================


#undef CNUM

#undef L_TAIL

#undef R_TAIL

#undef FLIP_SAMPLE

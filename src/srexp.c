#include "srexp.h"
#include "macro_var.h"

#define NAME srexp

#define CNUM 5

#define CNUM_SCALED 6

#define R_TAIL IT

#define R_ITF(u)(-(1.0 / (g.params[0])) * log(1 - (u)))

#define F(x)(g.params[0] * exp(-x * g.params[0]))

#define CDF(x)(1 - exp(-x * g.params[0]))

#define SCALABLE

#define SCALE(sample) sample / pp[0]

// =================================

#include "stors_sample.c"

// =================================

#include "stors_trunc_nav.c"

// =================================

#include "stors_sample_trunc.c"

// =================================


#undef CNUM

#undef R_TAIL

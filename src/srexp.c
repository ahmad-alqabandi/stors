#include "srexp.h"
#include "macro_var.h"

#define NAME srexp

#define CNUM 5

#define NON_SYMMETRIC_DIST

#define SCALE(sample) sample / pp[0]

#define R_TAIL IT

#define R_ITF(u)(-(1.0 / (g->params[0])) * log(1 - (u)))

#define F(x)(g->params[0] * exp(-x * g->params[0]))

#define CDF(x)(1 - exp(-x * g->params[0]))

// =================================

#define SCALABLE
#include "stors_sample_scalable_custom.c"
#undef SCALABLE
// =================================

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM
// =================================


#undef CNUM

#undef R_TAIL

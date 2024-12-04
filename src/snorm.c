#include "snorm.h"
#include "macro_var.h"

#define NAME srnorm

#define CNUM 1

#define SCALE(x) x * pp[1] + pp[0]

#define L_TAIL ARS

#define R_TAIL ARS

#define F(x) (1.0 / (g->params[1]  ) * exp(-0.5 * ((x - g->params[0]) / g->params[1]) * ((x - g->params[0]) / g->params[1])))

// =================================
# define SCALABLE
#include "stors_sample_scalable_custom.c"
#include "scaled_custom_check.c"
#undef SCALABLE

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#include "scaled_custom_check.c"
#undef CUSTOM
// =================================


#undef NAME

#define NAME srnorm_sym

#define FLIP_SAMPLE(sample, flip) flip ? g->symmetric - (sample - g->symmetric) : sample

# define SCALABLE
#include "stors_sample_scalable_custom.c"
#undef SCALABLE

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM

#undef FLIP_SAMPLE

#undef NAME

// =================================

#define NAME srnorm

#include "stors_trunc_nav.c"

#include "stors_sample_trunc.c"

// =================================

#undef CNUM

#undef L_TAIL

#undef R_TAIL
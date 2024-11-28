#include "snorm.h"
#include "macro_var.h"

#define NAME srnorm

#define CNUM 1

#define SCALE(sample) sample * pp[1] + pp[0]

#define L_TAIL ARS

#define R_TAIL ARS

#define F(sample) (1.0 / (g->params[1] * 2.50662827463) * exp(-0.5 * ((sample - g->params[0]) / g->params[1]) * ((sample - g->params[0]) / g->params[1])))

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

// #include "stors_sample_scalable.c"
// 
// #include "stors_sample_custom.c"

#undef FLIP_SAMPLE

#undef NAME

// =================================

#define NAME srnorm

#include "stors_trunc_nav.c"

#include "stors_sample_trunc.c"



#undef CNUM_SCALABLE

#undef CNUM_CUSTOM

#undef L_TAIL

#undef R_TAIL

#undef FLIP_SAMPLE

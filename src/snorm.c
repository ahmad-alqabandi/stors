#include "snorm.h"
#include "macro_var.h"

#define NAME srnorm

#define CNUM 1

#define SCALE(x) x * pp[1] + pp[0]

#define L_TAIL ARS

#define R_TAIL ARS

#define F(x) ( ( g->params[1]  ) * exp(-0.5 * ((x - g->params[0]) * g->params[1]) * ((x - g->params[0]) * g->params[1])))

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

#define INPLACE

# define SCALABLE
#include "stors_sample_scalable_custom.c"
#undef SCALABLE

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM

#undef INPLACE

// =================================

#define FLIP_SAMPLE(sample, flip) flip ? g->symmetric - (sample - g->symmetric) : sample


# define SCALABLE
#include "stors_sample_scalable_custom.c"
#undef SCALABLE

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM



#define INPLACE

#define FLIP_SAMPLE(sample, flip) flip ? g->symmetric - (sample - g->symmetric) : sample

# define SCALABLE
#include "stors_sample_scalable_custom.c"
#undef SCALABLE

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM

#undef INPLACE

#undef FLIP_SAMPLE

// =================================

#undef CNUM

#undef L_TAIL

#undef R_TAIL
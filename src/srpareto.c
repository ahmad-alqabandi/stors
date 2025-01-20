#include "srpareto.h"
#include "macro_var.h"

#define NAME srpareto

#define CNUM 13

#define NON_SYMMETRIC_DIST

#define R_TAIL IT

// params[0] Scale /x_m
// params[1] Shape /alpha
#define R_ITF(u)(g->params[0] * pow((1 - u),(-1 / g->params[1])))

#define F(x) ((g->params[1] * pow(g->params[0], g->params[1])) / pow(x, g->params[1] + 1))

#define CDF(x)(1 - pow((g->params[0] / x), g->params[1]))

// =================================

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#include "scaled_custom_check.c"
#undef CUSTOM

#define INPLACE

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM

#undef INPLACE


// =================================


#undef CNUM

#undef R_TAIL

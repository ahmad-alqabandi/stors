#include "srgamma.h"
#include "macro_var.h"
#include "math.h"

#define NAME srgamma

#define CNUM 6

#define L_TAIL ARS

#define R_TAIL ARS

#define SYMMETRIC FALSE

#define F(x) (1.0 / (tgamma(g.params[0]) * pow(g.params[1], g.params[0]))) * pow((x), (g.params[0]) - 1.0) * exp(-(x) / (g.params[1]))

// =================================

#include "stors_sample.c"

// =================================

#include "stors_trunc_nav.c"

// =================================


#include "stors_sample_trunc.c"


// =================================

#undef NAME

#undef CNUM

#undef L_TAIL

#undef R_TAIL

#undef SYMMETRIC

#undef F


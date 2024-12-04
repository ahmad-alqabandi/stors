#include "srchisq.h"
#include "macro_var.h"
#include "math.h"

#define NAME srchisq

#define CNUM 7

#define L_TAIL ARS

#define R_TAIL ARS

#define F(sample) ( (1.0 / ( pow(2.0, ((double)(g->params[0])/2.0)) * tgamma((double)(g->params[0])/2.0) )) \
* pow( (sample), ((double)(g->params[0])/2.0 - 1.0) ) * exp( - (sample)/2.0 ) )                                                \



// =================================

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM



// =================================

#undef NAME

#undef CNUM

#undef L_TAIL

#undef R_TAIL

#undef SYMMETRIC

#undef F


#include "snorm.h"
#include "macro_var.h"

#define NAME srnorm

#define CNUM 1

#define L_TAIL ARS

#define R_TAIL ARS

#define SYMMETRIC FALSE

// #define L_ITF(u)(log(2*u))

// #define R_ITF(u)(-log(2 - 2 * u))

//#define F(sample)(0.3989423 * exp(-0.5 * sample * sample))
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

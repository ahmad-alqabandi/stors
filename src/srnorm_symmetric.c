#include "srnorm_symmetric.h"
#include "macro_var.h"

#define NAME srnorm_symmetric

#define CNUM 4

#define L_TAIL ARS

#define R_TAIL ARS

#define SYMMETRIC TRUE

#define FLIP_SAMPLE(sample,flip)(flip ? - sample : sample)


// #define L_ITF(u)(log(2*u))

// #define R_ITF(u)(-log(2 - 2 * u))

#define F(sample)(0.3989423 * exp(-0.5 * sample * sample))

// =================================

#include "stors_sample.c"

// =================================



#undef CNUM

#undef L_TAIL

#undef R_TAIL

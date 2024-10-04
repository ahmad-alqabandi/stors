#include "srexp.h"
#include "macro_var.h"

#define NAME srexp

#define CNUM 3

#define R_TAIL IT

#define R_ITF(u)(-log(1 - u))

#define F(x)(exp(-x))

#define CDF(x)(1 - exp(-x))


// =================================

#include "stors_sample.c"

// =================================

#include "stors_trunc_nav.c"

// =================================

#include "stors_sample_trunc.c"

// =================================


#undef CNUM

#undef R_TAIL

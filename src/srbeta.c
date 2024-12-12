#include "srbeta.h"
#include "macro_var.h"
#include "math.h"

#define NAME srbeta

#define CNUM 11

#define NON_SYMMETRIC_DIST

#define L_TAIL ARS

#define R_TAIL ARS

#define F(sample) ( \
  (pow((sample), ((double)(g->params[0]) - 1.0)) * pow((1.0 - (sample)), ((double)(g->params[1]) - 1.0))) / \
  tgamma((double)(g->params[0])) / tgamma((double)(g->params[1])) * tgamma((double)(g->params[0]) + (double)(g->params[1])) \
)                                                              \
  
// =================================

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#include "scaled_custom_check.c"
#undef CUSTOM

// =================================

#define INPLACE

# define CUSTOM
#include "stors_sample_scalable_custom.c"
#undef CUSTOM

#undef INPLACE




// =================================

#undef NAME

#undef CNUM

#undef L_TAIL

#undef R_TAIL

#undef SYMMETRIC

#undef F


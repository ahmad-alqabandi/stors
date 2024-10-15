#include "stors.h"
#include "R_cache.h"
#include "snorm.h"
#include "srexp.h"
#include "srnorm_symmetric.h"
#include "laplace.h"
#include "old_srnorm.h"
#include "srchisq.h"
#include "srgamma.h"


static const R_CallMethodDef callMethods[]  = {
  {"stors", (DL_FUNC) &stors, -1},
  {"stors_trunc_nav", (DL_FUNC) &stors_trunc_nav, -1},
  {"stors_trunc", (DL_FUNC) &stors_trunc, -1},
  {"cache_grid", (DL_FUNC) &cache_grid, -1},
  {"free_cache", (DL_FUNC) &free_cache, -1},
  {"free_cache_cnum", (DL_FUNC) &free_cache_cnum, -1},
  {"srnorm", (DL_FUNC) &srnorm, -1},
  {"srnorm_trunc_nav", (DL_FUNC) &srnorm_trunc_nav, -1},
  {"srnorm_trunc", (DL_FUNC) &srnorm_trunc, -1},
  {"srnorm_symmetric", (DL_FUNC) &srnorm_symmetric, -1},
  {"old_srnorm", (DL_FUNC) &old_srnorm, -1},
  {"old_srnorm_trunc_nav", (DL_FUNC) &old_srnorm_trunc_nav, -1},
  {"old_srnorm_trunc", (DL_FUNC) &old_srnorm_trunc, -1},
  {"srexp", (DL_FUNC) &srexp, -1},
  {"srexp_trunc_nav", (DL_FUNC) &srexp_trunc_nav, -1},
  {"srexp_trunc", (DL_FUNC) &srexp_trunc, -1},
  {"laplace", (DL_FUNC) &laplace, -1},
  {"laplace_trunc_nav", (DL_FUNC) &laplace_trunc_nav, -1},
  {"laplace_trunc", (DL_FUNC) &laplace_trunc, -1},
  {"srchisq", (DL_FUNC) &srchisq, -1},
  {"srchisq_trunc_nav", (DL_FUNC) &srchisq_trunc_nav, -1},
  {"srchisq_trunc", (DL_FUNC) &srchisq_trunc, -1},
  {"srgamma", (DL_FUNC) &srgamma, -1},
  {"srgamma_trunc_nav", (DL_FUNC) &srgamma_trunc_nav, -1},
  {"srgamma_trunc", (DL_FUNC) &srgamma_trunc, -1},
  {NULL, NULL, 0}
};

void R_init_stors(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}




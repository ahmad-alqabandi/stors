#include "stors.h"
#include "R_cache.h"
#include "snorm.h"
#include "srexp.h"
#include "laplace.h"

static const R_CallMethodDef callMethods[]  = {
  {"stors", (DL_FUNC) &stors, -1},
  {"stors_trunc_nav", (DL_FUNC) &stors_trunc_nav, -1},
  {"stors_trunc", (DL_FUNC) &stors_trunc, -1},
  {"cache_grid", (DL_FUNC) &cache_grid, -1},
  {"free_cache", (DL_FUNC) &free_cache, -1},
  {"free_cache_cnum", (DL_FUNC) &free_cache_cnum, -1},
  {"print_cached_grids", (DL_FUNC) &print_cached_grids, -1},
  {"srnorm", (DL_FUNC) &srnorm, -1},
  {"srnorm_trunc_nav", (DL_FUNC) &srnorm_trunc_nav, -1},
  {"srnorm_trunc", (DL_FUNC) &srnorm_trunc, -1},
  {"srexp", (DL_FUNC) &srexp, -1},
  {"srexp_trunc_nav", (DL_FUNC) &srexp_trunc_nav, -1},
  {"srexp_trunc", (DL_FUNC) &srexp_trunc, -1},
  {"laplace", (DL_FUNC) &laplace, -1},
  {"laplace_trunc_nav", (DL_FUNC) &laplace_trunc_nav, -1},
  {"laplace_trunc", (DL_FUNC) &laplace_trunc, -1},
  {NULL, NULL, 0}
};

void R_init_stors(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}




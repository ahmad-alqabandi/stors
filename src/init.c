#include "stors.h"
#include "R_cache.h"
#include "snorm.h"

static const R_CallMethodDef callMethods[]  = {
  {"stors", (DL_FUNC) &stors, -1},
  {"cache_grid", (DL_FUNC) &cache_grid, -1},
  {"free_cache", (DL_FUNC) &free_cache, -1},
  {"print_cached_grids", (DL_FUNC) &print_cached_grids, -1},
  {"srnorm", (DL_FUNC) &srnorm, -1},
  {NULL, NULL, 0}
};


void R_init_stors(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}




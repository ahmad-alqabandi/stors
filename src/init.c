#include "stors.h"
#include "R_cache.h"
#include "snorm.h"
#include "srexp.h"
#include "laplace.h"
// #include "old_srnorm.h"
#include "srchisq.h"
#include "srgamma.h"
#include "grid_utils.h"
#include "err.h"


static const R_CallMethodDef callMethods[]  = {
  {"stors", (DL_FUNC) &stors, 4},
  {"stors_trunc_nav", (DL_FUNC) &stors_trunc_nav, 2},
  {"stors_trunc", (DL_FUNC) &stors_trunc, 3},
  {"grid_info", (DL_FUNC) &grid_info, 1},
  {"cache_grid", (DL_FUNC) &cache_grid, 15},
  {"free_cache", (DL_FUNC) &free_cache, 0},
  {"free_cache_cnum", (DL_FUNC) &free_cache_cnum, 1},
  {"grid_error", (DL_FUNC) &grid_error, 2},
  
  {"srnorm_scaled", (DL_FUNC) &srnorm_scaled, 2},
  {"srnorm_sym_scaled", (DL_FUNC) &srnorm_sym_scaled, 2},
  {"srnorm_scaled_check", (DL_FUNC) &srnorm_scaled_check, 2},
  
  {"srnorm_custom", (DL_FUNC) &srnorm_custom, 1},
  {"srnorm_sym_custom", (DL_FUNC) &srnorm_sym_custom, 1},
  {"srnorm_custom_check", (DL_FUNC) &srnorm_custom_check, 1},
  {"srnorm_trunc_nav", (DL_FUNC) &srnorm_trunc_nav, 3},
  {"srnorm_trunc", (DL_FUNC) &srnorm_trunc, 8},
  
  // {"old_srnorm", (DL_FUNC) &old_srnorm, 1},
  // {"old_srnorm_trunc_nav", (DL_FUNC) &old_srnorm_trunc_nav, 2},
  // {"old_srnorm_trunc", (DL_FUNC) &old_srnorm_trunc, 8},
  // {"srexp", (DL_FUNC) &srexp, 2},
  // {"srexp_trunc_nav", (DL_FUNC) &srexp_trunc_nav, 2},
  // {"srexp_trunc", (DL_FUNC) &srexp_trunc, 8},
  // {"laplace", (DL_FUNC) &laplace, 2},
  // {"laplace_trunc_nav", (DL_FUNC) &laplace_trunc_nav,2},
  // {"laplace_trunc", (DL_FUNC) &laplace_trunc, 8},
  // {"srchisq", (DL_FUNC) &srchisq, 2},
  // {"srchisq_trunc_nav", (DL_FUNC) &srchisq_trunc_nav, 2},
  // {"srchisq_trunc", (DL_FUNC) &srchisq_trunc, 7},
  // {"srgamma", (DL_FUNC) &srgamma, 2},
  // {"srgamma_trunc_nav", (DL_FUNC) &srgamma_trunc_nav, 2},
  // {"srgamma_trunc", (DL_FUNC) &srgamma_trunc, 7},
  
  {NULL, NULL, 0}
};

void R_init_stors(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}




#include "stors.h"
#include "R_cache.h"
#include "snorm.h"
#include "srexp.h"
#include "srlaplace.h"
// #include "old_srnorm.h"
#include "srchisq.h"
#include "srgamma.h"
#include "srpareto.h"
#include "grid_utils.h"
#include "srbeta.h"
#include "err.h"


static const R_CallMethodDef callMethods[]  = {
  {"stors", (DL_FUNC) &stors, 4},
  {"stors_trunc_nav", (DL_FUNC) &stors_trunc_nav, 2},
  {"stors_trunc", (DL_FUNC) &stors_trunc, 3},
  {"grid_info", (DL_FUNC) &grid_info, 1},
  {"cache_grid", (DL_FUNC) &cache_grid, 17},
  {"free_cache", (DL_FUNC) &free_cache, 0},
  {"free_cache_cnum", (DL_FUNC) &free_cache_cnum, 1},
  {"grid_error", (DL_FUNC) &grid_error, 2},


  // NORMAL
  {"srnorm_scaled", (DL_FUNC) &srnorm_scaled, 2},
  {"srnorm_scaled_inplace", (DL_FUNC) &srnorm_scaled_inplace, 2},

  {"srnorm_sym_scaled", (DL_FUNC) &srnorm_sym_scaled, 2},
  {"srnorm_sym_scaled_inplace", (DL_FUNC) &srnorm_sym_scaled_inplace, 2},
  {"srnorm_scaled_check", (DL_FUNC) &srnorm_scaled_check, 3},

  {"srnorm_custom", (DL_FUNC) &srnorm_custom, 1},
  {"srnorm_custom_inplace", (DL_FUNC) &srnorm_custom_inplace, 1},

  {"srnorm_sym_custom", (DL_FUNC) &srnorm_sym_custom, 1},
  {"srnorm_sym_custom_inplace", (DL_FUNC) &srnorm_sym_custom_inplace, 1},

  {"srnorm_custom_check", (DL_FUNC) &srnorm_custom_check, 2},

  // {"srnorm_trunc_nav", (DL_FUNC) &srnorm_trunc_nav, 3},
  // {"srnorm_trunc", (DL_FUNC) &srnorm_trunc, 8},


  // // LAPLACE
  {"srlaplace_scaled", (DL_FUNC) &srlaplace_scaled, 2},
  {"srlaplace_scaled_inplace", (DL_FUNC) &srlaplace_scaled_inplace, 2},

  {"srlaplace_sym_scaled", (DL_FUNC) &srlaplace_sym_scaled, 2},
  {"srlaplace_sym_scaled_inplace", (DL_FUNC) &srlaplace_sym_scaled_inplace, 2},
  {"srlaplace_scaled_check", (DL_FUNC) &srlaplace_scaled_check, 3},

  {"srlaplace_custom", (DL_FUNC) &srlaplace_custom, 1},
  {"srlaplace_custom_inplace", (DL_FUNC) &srlaplace_custom_inplace, 1},

  {"srlaplace_sym_custom", (DL_FUNC) &srlaplace_sym_custom, 1},
  {"srlaplace_sym_custom_inplace", (DL_FUNC) &srlaplace_sym_custom_inplace, 1},

  {"srlaplace_custom_check", (DL_FUNC) &srlaplace_custom_check, 2},

  // // EXPONENTIAL

  {"srexp_scaled", (DL_FUNC) &srexp_scaled, 2},
  {"srexp_scaled_inplace", (DL_FUNC) &srexp_scaled_inplace, 2},
  {"srexp_scaled_check", (DL_FUNC) &srexp_scaled_check, 3},

  {"srexp_custom", (DL_FUNC) &srexp_custom, 1},
  {"srexp_custom_inplace", (DL_FUNC) &srexp_custom_inplace, 1},
  {"srexp_custom_check", (DL_FUNC) &srexp_custom_check, 2},

  // // CHISQ

  {"srchisq_custom", (DL_FUNC) &srchisq_custom, 1},
  {"srchisq_custom_inplace", (DL_FUNC) &srchisq_custom_inplace, 1},
  {"srchisq_custom_check", (DL_FUNC) &srchisq_custom_check, 2},

  // // GAMMA

  {"srgamma_custom", (DL_FUNC) &srgamma_custom, 1},
  {"srgamma_custom_inplace", (DL_FUNC) &srgamma_custom_inplace, 1},
  {"srgamma_custom_check", (DL_FUNC) &srgamma_custom_check, 2},


  // BETA


  {"srbeta_custom", (DL_FUNC) &srbeta_custom, 1},
  {"srbeta_custom_inplace", (DL_FUNC) &srbeta_custom_inplace, 1},
  {"srbeta_custom_check", (DL_FUNC) &srbeta_custom_check, 2},


  // PARETO

  {"srpareto_custom", (DL_FUNC) &srpareto_custom, 1},
  {"srpareto_custom_inplace", (DL_FUNC) &srpareto_custom_inplace, 1},
  {"srpareto_custom_check", (DL_FUNC) &srpareto_custom_check, 2},

  {NULL, NULL, 0}
};

void R_init_stors(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}




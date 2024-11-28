
#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "cache.h"

#if defined(SCALABLE) || defined(CUSTOM) && defined(NAME) &&  defined(CNUM)

#ifdef SCALABLE
SEXP DEN_CHECK_SCALED(NAME)(SEXP s_size, SEXP Rpassed_params){
  struct grid *g = grids.grid + CNUM;
  
#endif
  
#ifdef CUSTOM
  SEXP DEN_CHECK_CUSTOM(NAME)(SEXP s_size){
    struct grid *g = grids.grid + CNUM + 1;
    
#endif
 
  
  if(g->x == NULL){
    REprintf("you need to optimize your destribution's grid first");
    R_RETURN_NULL;
  }
  
  if(g->is_symmetric == TRUE){
#ifdef SCALABLE
    return(DEN_SAMPLE_SYM_SCALED(NAME)(s_size,  Rpassed_params));
#endif
    
#ifdef CUSTOM
    return(DEN_SAMPLE_SYM_CUSTOM(NAME)(s_size));
#endif
    
  }else{
#ifdef SCALABLE
    return(DEN_SAMPLE_SCALED(NAME)(s_size, Rpassed_params));
#endif
    
#ifdef CUSTOM
    return(DEN_SAMPLE_CUSTOM(NAME)(s_size));
#endif

  } 
  
}
  
#endif

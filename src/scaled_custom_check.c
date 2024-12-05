
#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "cache.h"

#if defined(SCALABLE) || defined(CUSTOM) && defined(NAME) &&  defined(CNUM)

#ifdef SCALABLE
SEXP DEN_CHECK_SCALED(NAME)(SEXP s_size, SEXP Rpassed_params, SEXP R_reserved_memory){
  struct grid *g = grids.grid + CNUM;
  
#endif
  
#ifdef CUSTOM
  SEXP DEN_CHECK_CUSTOM(NAME)(SEXP s_size, SEXP R_reserved_memory){
    struct grid *g = grids.grid + CNUM + 1;
    
#endif
    
  
  if(g->x == NULL){
    REprintf("you need to optimize your destribution's grid first");
    R_RETURN_NULL;
  }
  
  if(g->is_symmetric == TRUE){
#ifdef SCALABLE
    if( R_reserved_memory == R_NilValue){
      return(DEN_SAMPLE_SYM_SCALED(NAME)(s_size,  Rpassed_params));
      
    }else{
      return(DEN_SAMPLE_SYM_SCALED_INPLACE(NAME)( Rpassed_params, R_reserved_memory));
      
    }
#endif
    
#ifdef CUSTOM
    if( R_reserved_memory == R_NilValue){
    return(DEN_SAMPLE_SYM_CUSTOM(NAME)(s_size));
      
    } else{
      return(DEN_SAMPLE_SYM_CUSTOM_INPLACE(NAME)(R_reserved_memory));
    }
#endif
    
  }else{
    
#ifdef SCALABLE
    if( R_reserved_memory == R_NilValue){
      return(DEN_SAMPLE_SCALED(NAME)(s_size, Rpassed_params));
    }else{
      return(DEN_SAMPLE_SCALED_INPLACE(NAME)(Rpassed_params, R_reserved_memory));
      
    }
#endif
    
#ifdef CUSTOM
    if( R_reserved_memory == R_NilValue){
    return(DEN_SAMPLE_CUSTOM(NAME)(s_size));
    } else{
      return(DEN_SAMPLE_CUSTOM_INPLACE(NAME)(R_reserved_memory));
    }
#endif

  }
  
}
  
#endif

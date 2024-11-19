#include "grid_utils.h"
#include "cache.h"

SEXP grid_info( SEXP Rcnum){
  
  int cnum = asInteger(Rcnum);
  
  if(grids.grid[cnum].x == NULL){
    R_RETURN_NULL
  }
  
  struct grid g = grids.grid[cnum];
  
  int offset = 1;
  
  int n = offset + g.n_params;
  
  SEXP grid_info = PROTECT(allocVector(REALSXP, n));
  
  double *grid_info_ptr = REAL(grid_info);
  
  grid_info_ptr[0] = g.is_symmetric;

  for(int i = offset; i < n; i++){
    
    grid_info_ptr[i] = g.params[i];
    
  }
  
  UNPROTECT(1);
  
  return grid_info;
  
}
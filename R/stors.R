#' Sampling Function for Users' Grid
#'
#' @description
#' This function generates a sampling function based on a grid created and optimized by the user using the `build_grid()` function. The resulting sampling function can then be used to produce samples.
#'
#' @param grid The sampling grid created by the user.
#'
#' @return 
#' Returns a sampling function that can be used to generate samples from the input `grid`.
#' 
#' @import digest digest
#' @export
stors <- function(grid) {
  
  force(grid)
  
  is_valid_grid(grid)
  
  if(digest(grid) %in% stors_env$user_cached_grids$Id ){
    
        Cnum = subset(stors_env$user_cached_grids, Id == digest(grid))$Cnum
    
  } else{
    
    n = nrow(stors_env$user_cached_grids) + 1
    
    stors_env$user_cached_grids[n,]$Id = digest(grid)
    
    Cnum = stors_env$user_cached_grids[n,]$Cnum = stors_env$grids$biultin$builtin_num + n - 1
    
    cache_grid_c(Cnum, grid)
    
  }
  
  dens_func <- eval(parse(text = grid$dens_func))
  
  rfunc_env <- new.env()
  
  rm(grid)
  
  function_string <- paste0("function(n) { .Call(C_stors,n,",paste0(Cnum),",dens_func,rfunc_env) }" )
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
}





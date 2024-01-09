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
  
  Cnum <- cache_stors_grid(grid)
  
  dens_func <- eval(parse(text = grid$dens_func))
  
  rfunc_env <- new.env()
  
  rm(grid)
  
  function_string <- paste0("function(n) { .Call(C_stors, n, ",paste0(Cnum),", dens_func,rfunc_env) }" )
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
}


#' @import digest digest
#' @export
trunc <- function(grid, xl, xr) {
  
  
  force(grid)
  
  is_valid_grid(grid)
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >=  grid$grid_bounds[1],
    "xr must be smaller than the density upper bound" = xr <=  grid$grid_bounds[2]
  )
  
  Cnum <- cache_stors_grid(grid)
  
  dens_func <- eval(parse(text = grid$dens_func))
  
  rfunc_env <- new.env()
  
  rm(grid)

  Upper_cumsum = .Call(C_stors_trunc_nav,Cnum, xl, xr)
  
  print(Upper_cumsum)
  
  function_string <- paste0("function(n) { .Call(C_stors_trunc, n, ",paste0(Cnum),", ",paste0(xl),", ",paste0(xr),", ",paste0(Upper_cumsum[1]),",", paste0(Upper_cumsum[2]),", ", paste0(as.integer(Upper_cumsum[3])),", ", paste0(as.integer(Upper_cumsum[4]))," , dens_func, rfunc_env) }")
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
  
  # return(
  #   function(n){
  #     .Call(C_stors_trunc, n ,Cnum, xl, xr, Upper_cumsum[1], Upper_cumsum[2], as.integer(Upper_cumsum[3]), as.integer(Upper_cumsum[4]), dens_func, rfunc_env)
  #   }
  # )

}


#' Title
#'
# d_stors_upper = function(n ,cnum, xl, xr ,csl , csr, il, ir, dens_func, rfunc_env){
#   
#   .Call(C_stors_trunc, n ,cnum, xl, xr, csl, csr, il, ir,  dens_func, rfunc_env)
#   
# }

cache_stors_grid <- function(grid){
  
  if(digest(grid) %in% stors_env$user_cached_grids$Id ){
    
    Cnum = subset(stors_env$user_cached_grids, Id == digest(grid))$Cnum
    
  } else{
    
    n = nrow(stors_env$user_cached_grids) + 5
    # 5 is padding
    
    stors_env$user_cached_grids[n,]$Id = digest(grid)
    
    Cnum <- stors_env$user_cached_grids[n,]$Cnum <- stors_env$grids$biultin$builtin_num + n
    
    cache_grid_c(Cnum, grid)
    
  }
  
  return(Cnum)
  
}

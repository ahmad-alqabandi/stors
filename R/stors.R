#' stros for sampling
#'
#' @param n sample size
#' @param grid proposal grid
#' @param f target density
#' @return sample of size n
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
    
    cash_grid_c(Cnum, grid)
    
  }
  
  dens_func <- eval(parse(text = grid$dens_func))
  
  rfunc_env <- new.env()
  
  rm(grid)
  
  function(n) {
    .Call(
      C_stors,
      n,
      eval(expression(Cnum)),
      dens_func,
      rfunc_env
    )
  }
}


#' stros for sampling
#'
#' @param n sample size
#' @param grid proposal grid
#' @param f target density
#' @return sample of size n
#' @export
#'
srnorm <- function(n, mean = 0, sd = 1) {
  # stopifnot("you need to optimize the grid first" = stors_env$grids_config$snorm$opt)
  .Call(C_srnorm, n)
}

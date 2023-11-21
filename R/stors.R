#' stros for sampling
#'
#' @param n sample size
#' @param grid proposal grid
#' @param f target density
#' @return sample of size n
#' @export
#' @import digest digest
stors <- function(grid) {
  force(grid)
  
  stopifnot(" This grid is not optmized using set_grid() " = digest(grid) %in% grids_env$grids_config$creatd_Id )
  
  n = nrow(grids_env$grids_config$grids)+1
  
  grids_env$grids_config$grids[n,]$Id = digest(grid)
  
  Cnum = grids_env$grids_config$grids[n,]$Cnum = grids_env$grids_config$builtin_num + n - 1
  
  grid$dens_func <- eval(parse(text = grid$dens_func))

  rfunc_env <- new.env()
  
  cash_grid_c(Cnum, grid)

  # MAKE SURE TO RM() ALL DATA AFTER IS GET CACHED IN C
  function(n) {
    .Call(
      C_stors,
      n,
      eval(expression(Cnum)),
      grid$dens_func,
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
snorm <- function(n, mean = 0, sd = 1) {
  # stopifnot("you need to optimize the grid first" = grids_env$grids_config$snorm$opt)
  .Call(C_snorm, n)
}

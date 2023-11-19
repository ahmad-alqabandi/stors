#' stros for sampling
#'
#' @param n sample size
#' @param grid proposal grid
#' @param f target density
#' @return sample of size n
#' @export
#'
stors <- function(dist_name) {
  force(dist_name)

  data_dir <- tools::R_user_dir("stors", "data")

  stopifnot(`there are no grids in user_data directory` = dir.exists(data_dir))

  grid_dir <- file.path(data_dir, paste0(dist_name, ".rds"))

  stopifnot(`there is no grid for 'dist_name' density` = file.exists(grid_dir))

  grid <- readRDS(grid_dir)

  grid$dens_func <- eval(parse(text = grid$dens_func))

  rfunc_env <- new.env()

  # MAKE SURE TO RM() ALL DATA AFTER IS GET CACHED IN C
  function(n) {
    .Call(
      C_stors,
      n,
      grid$grid_data$x,
      grid$grid_data$s_upper_lower,
      grid$grid_data$p_a,
      grid$steps_number,
      grid$sampling_probabilities,
      grid$unif_scaler,
      grid$grid_data$s_upper,
      grid$lt_properties,
      grid$rt_properties,
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
  .Call(C_snorm, n, 0)
}

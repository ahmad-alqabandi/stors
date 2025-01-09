#' Sampling Function for Users' Grid
#'
#' @description
#' This function generates a sampling function based on a grid created by the user using the \code{build_grid()} function.
#' The resulting sampling function can then be used to produce samples.
#'
#' @param grid The sampling grid created using the \code{build_grid()} function.
#'
#' @return
#' Returns a function that can be used to generate samples from the specified \code{grid}. If \code{xl} and \code{xr} are provided,
#' the samples are drawn from the truncated distribution within these bounds.
#'
#' @details
#' After a user creates a proposal grid for their desired sampling function using \code{\link{build_grid}},
#' this grid must be passed to \code{stors()} to create a sampling function for the target distribution.
#' \code{stors()} first checks whether the grid was indeed created using \code{build_grid()}. If the user has altered
#' or modified the grid returned from \code{build_grid()}, \code{stors()} will reject the altered grid; therefore,
#' no changes should be made to the grid after its creation. Once the grid is accepted by \code{stors()}, it is
#' cached in memory, allowing fast access to grid data for the compiled C code and reducing memory access latency.
#' Subsequently, \code{stors()} returns a function that can be utilized to generate samples from the target distribution,
#'
#' @examples
#'
#' # Example 1
#' # To sample from a standard normal distribution \( f(x) \sim \mathcal{N}(0,1) \),
#' # first build the proposal grid using \code{build_grid()}
#'
#' modes_norm = 0
#' f_norm <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x^2) }
#' h_norm <- function(x) { log(f_norm(x)) }
#' h_prime_norm <- function(x) { -x }
#' normal_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_norm,
#'  f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)
#'
#' # Generate samples from the standard normal distribution
#' sample_normal <- stors(normal_grid)
#' hist(sample_normal(100), main = "Normal Distribution Samples")
#'
#'
#' # Example 2
#' # Let's consider a bimodal distribution composed of two normal distributions:
#' # The first normal distribution N(0,1) with weight p = 0.3,
#' # and the second normal distribution N(4,1) with weight q = 0.7.
#'
#' f_bimodal <- function(x) {
#'  0.3 * (1 / sqrt(2 * pi) * exp(-0.5 * (x - 0)^2)) +
#'  0.7 * (1 / sqrt(2 * pi) * exp(-0.5 * (x - 4)^2))
#'}
#'
#' # Define the modes of the bimodal distribution
#'    modes_bimodal <- c(0.00316841, 3.99942)
#'
#' # Build the proposal grid for the bimodal distribution
#' bimodal_grid = build_grid(f = f_bimodal, modes = modes_bimodal, lb = -Inf, rb = Inf, steps = 1000)
#'
#' # Create the sampling function using \code{stors()}
#' sample_bimodal <- stors(bimodal_grid)
#'
#' # Generate and plot samples from the bimodal distribution
#' bimodal_samples <- sample_bimodal(1000)
#' hist(bimodal_samples, breaks = 30, main = "Bimodal Distribution Samples")
#'
#' # Create the truncated sampling function using \code{stors()} with truncation bounds [-0.5, 6]
#' truncated_bimodal_grid <- build_grid(f = f_bimodal, modes = modes_bimodal, lb = -0.5, rb = 6, steps = 1000)
#'
#' # Create the sampling function using \code{stors()}
#' sample_truncated_bimodal <- stors(truncated_bimodal_grid)
#'
#' # Generate and plot samples from the truncated bimodal distribution
#' truncated_sample <- sample_truncated_bimodal(1000)
#' hist(truncated_sample, breaks = 30, main = "Truncated Bimodal Distribution Samples")
#'
#'
#' @import digest digest
#' @export
stors <- function(grid) {

  force(grid)
  is_valid_grid(grid)
  c_num <- cache_user_grid_c(grid)

  # print(paste0("cnum = ", c_num))

  # if (xl != grid$grid_bounds[1] || xr != grid$grid_bounds[2]) {
  #   stopifnot(
  #     "xl must be a scaler" = (is.numeric(xl) && length(xl) == 1),
  #     "xr must be a scaler" = (is.numeric(xr) && length(xr) == 1),
  #     "xl must be smaller that xr" = xl < xr,
  #     "xl must be greater than or equal the density lower bound" = xl >  grid$grid_bounds[1],
  #     "xr must be smaller than or equal the density upper bound" = xr <  grid$grid_bounds[2]
  #   )
  #   upper_cumsum <- .Call(C_stors_trunc_nav, c_num, xl, xr)
  #   stopifnot(
  #     "xl is has a CDF close to 1" = (upper_cumsum[1] != 1),
  #     "xr is has a CDF close to 0" = (upper_cumsum[2] != 0)
  #   )
  #   function_string <- paste0("function(n) { .Call(C_stors_trunc, n, ", paste0(c_num), ", ", paste0(xl), ", ", paste0(xr), ", ", paste0(upper_cumsum[1]),
  #                             ",", paste0(upper_cumsum[2]), ", ", paste0(as.integer(upper_cumsum[3])), ", ", paste0(as.integer(upper_cumsum[4])),
  #                             " , dens_func, rfunc_env) }")
  # } else {
  #   function_string <- paste0("function(n) { .Call(C_stors, n, ", paste0(c_num), ", dens_func, rfunc_env) }")
  # }

  sampling_function <- parse(text = grid$dens_func)
  dens_func <- eval(sampling_function)

  rfunc_env <- new.env()

  function_string <- paste0("function(n) { .Call(C_stors, n, ", paste0(c_num), ", dens_func, rfunc_env) }")

  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)

  rm(grid)

  return(sampling_function)

}

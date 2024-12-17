# TODO finish all documentations. and make sure you are passing the check
# TODO get structure for the paper (look at computational paper structures) is jss
# TODO re-read the old STORS related papers
# TODO when you finish the documentation let Louis know via email so he can test it

# TODO watch Memento movie
# TODO watch Yesterday (2019) movie
# TODO watch dapperton bear movie
# TODO watch The Usual Suspects movie


#' Sampling from the Normal Distribution
#' @rdname srnorm
#' @order 1
#'
#' @description
#' The \code{srnorm()} function generates random samples from a Normal distribution using the STORS algorithm.
#' It employs an optimized proposal distribution around the mode and Adaptive Rejection Sampling (ARS) for the tails.
#'
#' @details
#'
#' The Normal distribution has the density:
#' \deqn{ f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mu)^2}{2\sigma^2}} }
#' where \eqn{\mu} is the mean and \eqn{\sigma} is the standard deviation.
#'
#' These two function are for sampling using the STORS algorithm based on the grid that has been constructed using \code{\link{srnorm_optimize}}.
#'
#' By default, \code{srnorm()} samples from a standard Normal distribution (\code{mean = 0}, \code{sd = 1}).
#' The proposal distribution is pre-optimized at package load time using \code{srnorm_optimize()} with
#' \code{steps = 4091}, creating a scalable grid centered around the mode.
#'
#' If \code{srnorm()} is called with custom \code{mean} or \code{sd} parameters, the samples are generated
#' from the standard Normal distribution, then scaled and location shifted accordingly.
#'
#' \section{When \code{srnorm_optimize()} is explicitly called:}{
#'\itemize{
#'  \item A grid is created and cached. If no parameters are provided, a standard grid is created (\code{mean = 0}, \code{sd = 1}).
#'  \item Providing \code{mean} or \code{sd} creates a custom grid, which is cached for use with \code{srnorm_custom()}.
#'  \item The optimization process can be controlled via parameters such as \code{steps}, \code{grid_range}, or
#'   \code{theta}. If no parameters are provided, the grid is optimized via brute force based on the.
#'   \code{target_sample_size}.
#'}
#'}
#'
#' For efficiency, \code{srnorm_optimize()} also supports one-tailed grids via the \code{symmetric} argument.
#' This reduces memory usage for symmetrical distributions, though sampling may be slower for large samples
#' due to the additional computation required for symmetry checks.
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param mean Numeric. Mean parameter of the Normal distribution.
#' @param sd Numeric. Standard deviation of the target Normal distribution.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is over
#' written in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing samples from the Normal distribution with the specified
#' #TODO this is ok but when splitting. \describe{ \item{. ...}} for the grid optimizer
#' \code{mean} and \code{sd}.
#'
#'**NOTE:** when the x parameter is specified then for performance reason x is updated in-place with the simulation
#'
#'@seealso
#'
#' \code{\link{srnorm_optimize}} to optimize the custom or the scaled proposal grid.
#'
#' @examples
#' # Generate 10 samples from the standard Normal distribution
#' samples <- srnorm(10)
#' print(samples)
#'
#' # Generate 10 samples using a pre-allocated vector
#' x <- numeric(10)
#' srnorm(10, x = x)
#' print(x)
#'
#' # Generate 10 samples from a Normal distribution with mean = 2 and sd = 3
#' samples <- srnorm(10, mean = 2, sd = 3)
#' print(samples)
#'
#' @export
srnorm <- function(n = 1, mean = 0, sd = 1, x = NULL) {
  .Call(C_srnorm_scaled_check, n, c(mean, sd), x)
}


#' Sampling from Custom Normal Distribution
#' @rdname srnorm
#' @order 2
#'
#' @export
srnorm_custom <- function(n = 1, x = NULL) {
  .Call(C_srnorm_custom_check, n, x)
}


#' Optimizing Normal Distribution Grid
#'
#' @param mean Numeric. Mean parameter of the Normal distribution.
#' @param sd Numeric. Standard deviation of the target Normal distribution.
#' @param xl (optional) Scalar lower bounds.
#' @param xr (optional) Scalar upper bounds.
#' @param steps (optional) Scalar integer indicating the number of steps in the proposal distribution.
#' @param grid_range (optional) Optional Vector of two elements specifying the start and end points for constructing steps along the x-axis.
#' @param theta (optional) the minimum acceptable pre-acceptance value for the grid which indirectly determine the number of steps. Can not be specified at the same time as `steps`.
#' @param target_sample_size (optional) Scalar integer indicating the target sample size optimizer will simulate batches of this size when benchmarking.
#' @param verbose (optional) Logical if set to `TRUE`, a table detailing the optimization areas and steps will be displayed during grid optimization.
#' @param symmetric (optional) Logical if True then only one half of the grid is built ti the left to the mode and a coin toss is used to determine weather to reflect around the mode, this allows a dencer grid for the same memory at the cost of an additional uniform draw "default is `FALSE`"
#'
#' @export
srnorm_optimize <- function(
  mean = NULL,
  sd = NULL,
  xl = -Inf,
  xr = Inf,
  steps = NULL,
  grid_range = NULL,
  theta = NULL,
  target_sample_size = 1000,
  verbose = FALSE,
  symmetric = NULL) {

  dist_name <- "srnorm"

  dendata <- pbgrids[[dist_name]]

  f_params <- list(mean = mean, sd = sd)

  if (dendata$scalable) {
    isnull <- sapply(f_params, is.null)

    if (all(isnull)) {
      cnum <- dendata$c_num
      grid_type <- "scaled"
    } else {
      cnum <- dendata$c_num + 1
      grid_type <- "custom"
    }

    f_params <- ifelse(isnull, dendata$std_params, f_params)

  } else {
    cnum <- dendata$c_num + 1
    grid_type <- "custom"
  }

  modes <- dendata$set_modes(f_params$mean)

  f <- dendata$create_f(f_params$mean, f_params$sd)


  check_grid_optimization_criteria(symmetric, cnum, dendata)

  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type, symmetric, cnum, verbose)

}

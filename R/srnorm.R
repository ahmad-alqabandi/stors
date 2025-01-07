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
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param mean Numeric. Mean parameter of the Normal distribution.
#' @param sd Numeric. Standard deviation of the target Normal distribution.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is over
#' written in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing samples from the Normal distribution with the specified
#' \code{mean} and \code{sd}.
#'
#' \bold{NOTE:} When the \code{x} parameter is specified, it is updated in-place with the simulation for performance reasons.
#'
#' @seealso
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
#' @description
#' The \code{srnorm_optimize()} function generates an optimized proposal grid for a targeted Normal distribution.
#'  The grid can be customized and adjusted based on various options provided by the user.
#'
#'
#' @details
#'When \code{srnorm_optimize()} is explicitly called:
#'\itemize{
#'  \item A grid is created and cached. If no parameters are provided, a standard grid is created (\code{mean = 0}, \code{sd = 1}).
#'  \item Providing \code{mean} or \code{sd} creates a custom grid, which is cached for use with \code{srnorm_custom()}.
#'  \item The optimization process can be controlled via parameters such as \code{steps}, \code{grid_range}, or
#'   \code{theta}. If no parameters are provided, the grid is optimized via brute force based on the.
#'   \code{target_sample_size}.
#'}
#'
#' @param mean (optional) Numeric. Mean parameter of the Normal distribution. Defaults to \code{NULL}, which implies a scalable grid with \code{mean = 0}.
#' @param sd (optional) Numeric. Standard deviation of the target Normal distribution. Defaults to \code{NULL}, which implies a scalable grid with \code{sd = 1}.
#' @param xl Numeric. Left truncation bound for the target distribution. Defaults to \code{-Inf}, representing no left truncation.
#' @param xr Numeric. Right truncation bound for the target distribution. Defaults to \code{Inf}, representing no right truncation.
#' @param steps (optional) Integer. Desired number of steps in the proposal grid. Defaults to \code{NULL}, which means the number of steps is determined automatically during optimization.
#' @param grid_range (optional) Numeric vector. Specifies the range for optimizing the steps part of the proposal grid. Defaults to \code{NULL}, indicating automatic range selection.
#' @param theta (optional) Numeric. A parameter for grid optimization. Defaults to \code{NULL}.
#' @param target_sample_size (optional) Integer. Target sample size for grid optimization. Defaults to \code{1000}.
#' @param verbose Boolean. If \code{TRUE}, detailed optimization information, including areas and steps, will be displayed. Defaults to \code{FALSE}.
#' @param symmetric Boolean. If \code{TRUE}, the proposal will target only the right tail of the distribution, reducing the size of the cached grid and making sampling more memory-efficient.
#'  An additional uniform random number will be sampled to determine the sample's position relative to the mode of the distribution.
#'   While this improves memory efficiency, the extra sampling may slightly impact performance, especially when the grid efficiency is close to 1. Defaults to \code{FALSE}.
#'
#'
#' @return
#' The user does not need to store the returned value, because the package internally cashes the grid. However, we explain here the full returned grid for advanced users.
#'
#' A list containing the optimized grid and related parameters for the specified built-in distribution:
#' \describe{
#'   \item{\code{grid_data}}{A data frame with detailed information about the grid steps, including:
#'   \describe{
#'     \item{\code{x}}{The start point of each step on the x-axis.}
#'     \item{\code{s_upper}}{The height of each step on the y-axis.}
#'     \item{\code{p_a}}{Pre-acceptance probability for each step.}
#'     \item{\code{s_upper_lower}}{A vector used to scale the uniform random number when the sample is accepted.}
#'   }}
#'   \item{\code{areas}}{A numeric vector containing the areas under:
#'   \describe{
#'     \item{\code{left_tail}}{The left tail bound.}
#'     \item{\code{steps}}{The middle steps.}
#'     \item{\code{right_tail}}{The right tail bound.}
#'   }}
#'   \item{\code{steps_number}}{An integer specifying the number of steps in the proposal.}
#'   \item{\code{sampling_probabilities}}{A numeric vector with:
#'   \describe{
#'     \item{\code{left_tail}}{The probability of sampling from the left tail.}
#'     \item{\code{left_and_middle}}{The combined probability of sampling from the left tail and middle steps.}
#'   }}
#'   \item{\code{unif_scaler}}{A numeric scalar, the inverse probability of sampling from the steps part of the proposal (\eqn{\frac{1}{p(lb < x < rb)}}). Used for scaling uniform random values.}
#'   \item{\code{lt_properties}}{A numeric vector of 5 values required for Adaptive Rejection Sampling (ARS) in the left tail.}
#'   \item{\code{rt_properties}}{A numeric vector of 6 values required for ARS in the right tail.}
#'   \item{\code{alpha}}{A numeric scalar representing the uniform step area.}
#'   \item{\code{tails_method}}{A string, either \code{"ARS"} (Adaptive Rejection Sampling) or \code{"IT"} (Inverse Transform), indicating the sampling method for the tails.}
#'   \item{\code{grid_bounds}}{A numeric vector specifying the left and right bounds of the target density.}
#'   \item{\code{cnum}}{An integer representing the cache number of the created grid in memory.}
#'   \item{\code{symmetric}}{A numeric scalar indicating the symmetry point of the grid, or \code{NULL} if not symmetric.}
#'   \item{\code{f_params}}{A list of parameters for the target density that the proposal grid is designed for.}
#'   \describe{
#'     \item{\code{mean}}{The mean of the target distribution.}
#'     \item{\code{sd}}{The standard deviation of the target distribution.}
#'   }
#'   \item{\code{is_symmetric}}{A logical value indicating whether the proposal grid is symmetric.}
#'   \item{\code{grid_type}}{A string indicating the type of the generated grid:
#'   \describe{
#'     \item{\code{"scaled"}}{The grid is "scalable" and standardized with \code{mean = 0} and \code{sd = 1}. This is used when parameters \code{mean} and \code{sd} are either \code{NULL} or not provided. Scalable grids are compatible with \code{\link{srnorm}}.}
#'     \item{\code{"custom"}}{The grid is "custom" when either \code{mean} or \code{sd} is provided. Custom grids are compatible with \code{\link{srnorm_custom}}.}
#'   }}
#'   \item{\code{target_function_area}}{A numeric scalar estimating the area of the target distribution.}
#'   \item{\code{dens_func}}{A string containing the hardcoded density function.}
#'   \item{\code{density_name}}{A string specifying the name of the target density distribution.}
#'   \item{\code{lock}}{An identifier used for saving and loading the grid from disk.}
#' }
#'
#' @seealso
#' \code{\link{srnorm}}: Function to sample from a scalable grid generated by \code{srnorm_optimize()}.
#' \code{\link{srnorm_custom}}: Function to sample from a custom grid tailored to user specifications.
#'
#'
#' @examples
#' # Generate scalable grid that with mean = 0 and sd = 1, that has 4096 steps
#' scalable_grid <- srnorm_optimize(steps = 4096)
#'
#' # Generate custom grid that with mean = 2 and sd = 1
#' scalable_grid <- srnorm_optimize(mean = 2, sd = 1)
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
  symmetric = FALSE) {

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

  if (symmetric) {
    symmetric <- modes
  } else {
    symmetric <- NULL
    }

  f <- dendata$create_f(f_params$mean, f_params$sd)


  check_grid_opt_criteria(symmetric, cnum, dendata)

  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type, symmetric, cnum, verbose)

}

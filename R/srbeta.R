
#' Sampling from the Beta Distribution
#' @rdname srbeta_custom
#'
#' @description
#' The \code{srbeta_custom()} function generates random samples from a Beta distribution using the STORS algorithm.
#' It employs an optimized proposal distribution around the mode and Adaptive Rejection Sampling (ARS) for the tails.
#'
#' @details
#'
#' The Beta Distribution
#' @details
#' The Beta distribution has the probability density function (PDF):
#' \deqn{f(x | \alpha, \beta) = \frac{\Gamma(\alpha + \beta)}{\Gamma(\alpha)\Gamma(\beta)} x^{\alpha - 1} (1 - x)^{\beta - 1}, \quad 0 \leq x \leq 1,}
#' where:
#' \itemize{
#'   \item{\eqn{\alpha}}{is the first shape parameter (\eqn{\alpha > 0}).}
#'   \item{\eqn{\beta}}{is the second shape parameter (\eqn{\beta > 0}).}
#' }
#' The Beta distribution is widely used in Bayesian statistics and in modeling probabilities and proportions.
#'
#' This function samples from a grid constructed using \code{\link{srbeta_optimize}}, employing the STORS algorithm.
#'
#' By default, \code{srbeta_custom()} samples from the standard Beta distribution with \code{shape1 = 1} and \code{shape2 = 1}.
#' The proposal distribution is pre-optimized at package load time using \code{srbeta_optimize()} with
#' \code{steps = 4091}, creating a scalable grid centered around the mode.
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param shape1 Numeric. The first shape parameter (\eqn{\alpha}) of the Beta distribution.
#' @param shape2 Numeric. The second shape parameter (\eqn{\beta}) of the Beta distribution.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is overwritten in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing random samples from the Beta distribution.
#' The \code{shape1} and \code{shape2} parameters are specified during the optimization process using \code{srbeta_optimize()}.
#'
#' \bold{NOTE:} When the \code{x} parameter is specified, it is updated in-place with the simulation for performance reasons.
#'
#' @seealso
#' \code{\link{srbeta_optimize}} to optimize the custom grid.
#'
#' @examples
#' # Generate 10 samples from Beta Distribution
#' samples <- srbeta_custom(10)
#' print(samples)
#'
#' # Generate 10 samples using a pre-allocated vector
#' x <- numeric(10)
#' srbeta_custom(10, x = x)
#' print(x)
#'
#' @export
srbeta_custom <- function(n = 1, x = NULL) {
  .Call(C_srbeta_custom_check, n, x)
}





#' Optimizing Beta Distribution Grid
#' @description
#' The \code{srbeta_optimize()} function generates an optimized proposal grid for a targeted Beta distribution.
#' The grid can be customized and adjusted based on various options provided by the user.
#'
#' @details
#' When \code{srbeta_optimize()} is explicitly called:
#' \itemize{
#'   \item A grid is created and cached. If no parameters are provided, a standard grid is created with \code{shape1 = 1} and \code{shape2 = 1}.
#'   \item Providing \code{shape1} and \code{shape2} creates a custom grid, which is cached for use with \code{srbeta_custom()}.
#'   \item The optimization process can be controlled via parameters such as \code{steps}, \code{grid_range}, or
#'   \code{theta}. If no parameters are provided, the grid is optimized via brute force based on the
#'   \code{target_sample_size}.
#' }
#'
#' @param shape1 (optional) Numeric. The first shape parameter (\eqn{\alpha}) of the Beta distribution. Defaults to \code{1}.
#' @param shape2 (optional) Numeric. The second shape parameter (\eqn{\beta}) of the Beta distribution. Defaults to \code{1}.
#' @param xl Numeric. Left truncation bound for the target distribution. Defaults to \code{0}, as the Beta distribution is defined only on the interval \code{[0, 1]}.
#' @param xr Numeric. Right truncation bound for the target distribution. Defaults to \code{1}, as the Beta distribution is defined only on the interval \code{[0, 1]}.
#' @param steps (optional) Integer. Desired number of steps in the proposal grid. Defaults to \code{4091}.
#' @param grid_range (optional) Numeric vector. Specifies the range for optimizing the steps part of the proposal grid. Defaults to \code{NULL}, indicating automatic range selection.
#' @param theta (optional) Numeric. A parameter for grid optimization. Defaults to \code{NULL}.
#' @param target_sample_size (optional) Integer. Target sample size for grid optimization. Defaults to \code{1000}.
#' @param verbose Boolean. If \code{TRUE}, detailed optimization information, including areas and steps, will be displayed. Defaults to \code{FALSE}.
#'
#' @return
#' A list containing the optimized grid and related parameters for the specified Beta distribution. The grid is also cached for internal use.
#' \describe{
#'   \item{\code{grid_data}}{Detailed information about the grid steps, including \code{x}, \code{s_upper}, \code{p_a}, and \code{s_upper_lower}.}
#'   \item{\code{areas}}{The areas under the left tail, steps, and right tail of the proposal distribution.}
#'   \item{\code{steps_number}}{The number of steps in the proposal grid.}
#'   \item{\code{f_params}}{The parameters (\code{shape1} and \code{shape2}) of the Beta distribution.}
#' }
#'
#' @seealso
#' \code{\link{srbeta_custom}}: Function to sample from a custom grid generated by \code{srbeta_optimize()}.
#'
#' @examples
#' # Generate a standard grid with shape1 = 1 and shape2 = 1
#' standard_grid <- srbeta_optimize()
#'
#' # Generate a custom grid with shape1 = 2 and shape2 = 3
#' custom_grid <- srbeta_optimize(shape1 = 2, shape2 = 3)
#'
#' @export
srbeta_optimize <- function(
    shape1 = NULL,
    shape2 = NULL,
    xl = 0,
    xr = 1,
    steps = NULL,
    grid_range = NULL,
    theta = NULL,
    target_sample_size = 1000,
    verbose = FALSE) {

  symmetric <- NULL

  dist_name <- "srbeta"

  dendata <- pbgrids[[dist_name]]

  if (shape1 <= 1 || shape2 <= 1) {
    message("Grid building is not available for shape1 <= 1 or shape2 <= 1.")
    return()
  }

  f_params <- list(shape1 = shape2, shape2 = shape2)

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

  modes <- dendata$set_modes(f_params$shape1, f_params$shape2)

  f <- dendata$create_f(f_params$shape1, f_params$shape2)

  check_grid_opt_criteria(symmetric, cnum, dendata)

  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type, symmetric, cnum, verbose)

}

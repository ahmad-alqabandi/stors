#' Sampling from Pareto Distribution
#' @rdname srpareto_custom
#'
#' @description
#' The \code{srpareto_custom()} function generates random samples from a Pareto distribution using the STORS algorithm.
#' It employs an optimized proposal distribution around the mode and Adaptive Rejection Sampling (ARS) for the tails.
#'
#' @details
#'
#' The Pareto Distribution
#' @details
#' The Pareto distribution has the probability density function (PDF):
#' \deqn{f(x | \alpha, \beta) = \frac{\beta^\alpha}{\Pareto(\alpha)} x^{\alpha - 1} \exp(-\beta x), \quad x \geq 0,}
#' where:
#' \describe{
#'   \item{\eqn{\alpha}}{is the shape parameter (\eqn{\alpha > 0}), which determines the shape of the distribution.}
#'   \item{\eqn{\beta}}{is the rate parameter (\eqn{\beta > 0}), which determines the rate of decay.}
#' }
#' The Pareto distribution is widely used in statistics, particularly in Bayesian inference and modeling waiting times.
#'
#' This function samples from a grid constructed using \code{\link{srpareto_optimize}}, employing the STORS algorithm.
#'
#' By default, \code{srpareto_custom()} samples from the standard Pareto distribution with \code{shape = 1} and \code{rate = 1}.
#' The proposal distribution is pre-optimized at package load time using \code{srpareto_optimize()} with
#' \code{steps = 4091}, creating a scalable grid centered around the mode.
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is overwritten in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing random samples from the Pareto distribution.
#' The \code{shape} and \code{rate} parameters are specified during the optimization process using \code{srpareto_optimize()}.
#'
#' \bold{NOTE:} When the \code{x} parameter is specified, it is updated in-place with the simulation for performance reasons.
#'
#' @note
#' This function is not scalable. Therefore, only the \code{srpareto_custom()} version is available, which requires the grid to be pre-optimized using \code{srpareto_optimize()} before calling this function.
#'
#' @seealso
#' \code{\link{srpareto_optimize}} to optimize the custom grid.
#'
#' @examples
#' # Generate 10 samples from Pareto Distribution
#' samples <- srpareto_custom(10)
#' print(samples)
#'
#' # Generate 10 samples using a pre-allocated vector
#' x <- numeric(10)
#' srpareto_custom(10, x = x)
#' print(x)
#'
#' @export
srpareto_custom <- function(n = 1, x = NULL) {
  .Call(C_srpareto_custom_check, n, x)
}



#' Optimizing Pareto Distribution Grid
#' @description
#' The \code{srpareto_optimize()} function generates an optimized proposal grid for a targeted Pareto Distribution.
#'  The grid can be customized and adjusted based on various options provided by the user.
#'
#'
#' @details
#'When \code{srpareto_optimize()} is explicitly called:
#'\itemize{
#'  \item A grid is created and cached. If no parameters are provided, a standard grid is created with \code{rate = 1}.
#'  \item Providing \code{rate} creates a custom grid, which is cached for use with \code{srpareto_custom()}.
#'  \item The optimization process can be controlled via parameters such as \code{steps}, \code{grid_range}, or
#'   \code{theta}. If no parameters are provided, the grid is optimized via brute force based on the.
#'   \code{target_sample_size}.
#'}
#'
#' @param rate (optional) Numeric. rate parameter of the Pareto Distribution. Defaults to \code{NULL}, which implies a scalable grid with \code{rate = 1}.
#' @param xl Numeric. Left truncation bound for the target distribution. Defaults to \code{-Inf}, representing no left truncation.
#' @param xr Numeric. Right truncation bound for the target distribution. Defaults to \code{Inf}, representing no right truncation.
#' @param steps (optional) Integer. Desired number of steps in the proposal grid. Defaults to \code{NULL}, which means the number of steps is determined automatically during optimization.
#' @param grid_range (optional) Numeric vector. Specifies the range for optimizing the steps part of the proposal grid. Defaults to \code{NULL}, indicating automatic range selection.
#' @param theta (optional) Numeric. A parameter for grid optimization. Defaults to \code{NULL}.
#' @param target_sample_size (optional) Integer. Target sample size for grid optimization. Defaults to \code{1000}.
#' @param verbose Boolean. If \code{TRUE}, detailed optimization information, including areas and steps, will be displayed. Defaults to \code{FALSE}.
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
#'     \item{\code{rate}}{ the rate of the target distribution.}
#'   }
#'   \item{\code{is_symmetric}}{A logical value indicating whether the proposal grid is symmetric.}
#'   \item{\code{grid_type}}{A string indicating the type of the generated grid:
#'   \describe{
#'     \item{\code{"scaled"}}{The grid is "scalable" and standardized with \code{rate = 1}. This is used when parameter \code{rate} is either \code{NULL} or not provided. Scalable grids are compatible with \code{\link{srpareto}}.}
#'     \item{\code{"custom"}}{The grid is "custom" when \code{rate} is provided. Custom grids are compatible with \code{\link{srpareto_custom}}.}
#'   }}
#'   \item{\code{target_function_area}}{A numeric scalar estimating the area of the target distribution.}
#'   \item{\code{dens_func}}{A string containing the hardcoded density function.}
#'   \item{\code{density_name}}{A string specifying the name of the target density distribution.}
#'   \item{\code{lock}}{An identifier used for saving and loading the grid from disk.}
#' }
#'
#' @seealso
#' \code{\link{srpareto}}: Function to sample from a scalable grid generated by \code{srpareto_optimize()}.
#' \code{\link{srpareto_custom}}: Function to sample from a custom grid tailored to user specifications.
#'
#'
#' @examples
#' # Generate scalable grid that with rate = 1, that has 4096 steps
#' scalable_grid <- srpareto_optimize(steps = 4096)
#'
#' # Generate custom grid that with rate = 4
#' scalable_grid <- srpareto_optimize(rate = 4)
#'
#' @export
srpareto_optimize <- function(scale = NULL,
                              shape = NULL,
                              xl = NULL,
                              xr = NULL,
                              steps = 4091,
                              grid_range = NULL,
                              theta = NULL,
                              target_sample_size = 1000,
                              verbose = FALSE) {
  dist_name <- "srpareto"

  dendata <- pbgrids[[dist_name]]

  f_params <- list(scale = scale, shape = shape)

  cnum <- dendata$c_num + 1

  grid_type <- "custom"

  f_params <- ifelse(sapply(f_params, is.null), dendata$std_params, f_params)

  modes <- dendata$set_modes(f_params$scale)

  symmetric <- NULL

  f <<- dendata$create_f(f_params$scale, f_params$shape)

  #TODO: exception for this distribution. this could produce issue.
  dendata$lb <- f_params$scale

  grid_optimizer(
    dendata,
    dist_name,
    xl,
    xr,
    f,
    modes,
    f_params,
    steps,
    grid_range,
    theta,
    target_sample_size,
    grid_type,
    symmetric,
    cnum,
    verbose
  )

}

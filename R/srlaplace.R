#' Sampling from Laplace Distribution
#' @rdname srlaplace
#' @order 1
#'
#' @description
#' The \code{srlaplace()} function generates random samples from a Laplace Distribution using the STORS algorithm.
#' It employs an optimized proposal distribution around the mode and Inverse Transform (IT) method for the tails.
#'
#' @details
#'
#' The Laplace distribution has the probability density function (PDF):
#' \eqn{f(x | \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right),}
#' where:
#' \describe{
#'   \item{\eqn{\mu}}{is the location parameter (mean of the distribution).}
#'   \item{\code{b}}{is the scale parameter, which controls the spread of the distribution (\code{b > 0}).}
#' }
#'
#' These two functions are for sampling using the STORS algorithm based on the proposal that has been constructed using \code{\link{srlaplace_optimize}}.
#'
#' By default, \code{srlaplace()} samples from a standard Laplace Distribution (\code{mu = 0}, \code{b = 1}).
#' The proposal distribution is pre-optimized at package load time using \code{srlaplace_optimize()} with
#' \code{steps = 4091}, creating a scalable proposal centered around the mode.
#'
#' If \code{srlaplace()} is called with custom \code{mu} or \code{b} parameters, the samples are generated
#' from the standard Laplace Distribution, then scaled and location shifted accordingly.
#'
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param mu Numeric, location parameter.
#' @param b Numeric, scale parameter.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is over
#' written in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing samples from the Laplace Distribution with the specified
#' \code{mu} and \code{b}.
#'
#' \bold{NOTE:} When the \code{x} parameter is specified, it is updated in-place with the simulation for performance reasons.
#'
#' @seealso
#' \code{\link{srlaplace_optimize}} to optimize the custom or the scaled proposal.
#'
#' @examples
#' # Generate 10 samples from the standard Laplace Distribution
#' samples <- srlaplace(10)
#' print(samples)
#'
#' # Generate 10 samples using a pre-allocated vector
#' x <- numeric(10)
#' srlaplace(10, x = x)
#' print(x)
#'
#' # Generate 10 samples from a Laplace Distribution with mu = 2 and b = 3
#' samples <- srlaplace(10, mu = 2, b = 3)
#' print(samples)
#'
#' @export
srlaplace <- function(n = 1, mu = 0, b = 1, x = NULL) {
  .Call(C_srlaplace_scaled_check, n, c(mu, b), x)
}


#' #' Sampling from Custom Laplace Distribution
#' @rdname srlaplace
#' @order 2
#'
#' @export
srlaplace_custom <- function(n = 1, x = NULL) {
  .Call(C_srlaplace_custom_check, n, x)
}



#' Optimizing Laplace Distribution proposal
#' @description
#' The \code{srlaplace_optimize()} function generates an optimized proposal for a targeted Laplace Distribution.
#'  The proposal can be customized and adjusted based on various options provided by the user.
#'
#'
#' @details
#'When \code{srlaplace_optimize()} is explicitly called:
#'\itemize{
#'  \item A proposal is created and cached. If no parameters are provided, a standard proposal is created (\code{mu = 0}, \code{b = 1}).
#'  \item Providing \code{mu} or \code{b} creates a custom proposal, which is cached for use with \code{srlaplace_custom()}.
#'  \item The optimization process can be controlled via parameters such as \code{steps}, \code{proposal_range}, or
#'   \code{theta}. If no parameters are provided, the proposal is optimized via brute force based on the.
#'   \code{target_sample_size}.
#'}
#'
#' @param mu (optional) Numeric, location parameter.
#' @param b (optional) Numeric, scale parameter.
#' @param xl Numeric. Left truncation bound for the target distribution. Defaults to \code{-Inf}, representing no left truncation.
#' @param xr Numeric. Right truncation bound for the target distribution. Defaults to \code{Inf}, representing no right truncation.
#' @param steps (optional) Integer. Desired number of steps in the proposal. Defaults to \code{NULL}, which means the number of steps is determined automatically during optimization.
#' @param proposal_range (optional) Numeric vector. Specifies the range for optimizing the steps part of the proposal. Defaults to \code{NULL}, indicating automatic range selection.
#' @param theta (optional) Numeric. A parameter for proposal optimization. Defaults to \code{NULL}.
#' @param target_sample_size (optional) Integer. Target sample size for proposal optimization. Defaults to \code{1000}.
#' @param verbose Boolean. If \code{TRUE}, detailed optimization information, including areas and steps, will be displayed. Defaults to \code{FALSE}.
#' @param symmetric Boolean. If \code{TRUE}, the proposal will target only the right tail of the distribution, reducing the size of the cached proposal and making sampling more memory-efficient.
#'  An additional uniform random number will be sampled to determine the sample's position relative to the mode of the distribution.
#'   While this improves memory efficiency, the extra sampling may slightly impact performance, especially when the proposal efficiency is close to 1. Defaults to \code{FALSE}.
#'
#'
#' @return
#' The user does not need to store the returned value, because the package internally cashes the proposal. However, we explain here the full returned proposal for advanced users.
#'
#' A list containing the optimized proposal and related parameters for the specified built-in distribution:
#' \describe{
#'   \item{\code{data}}{A data frame with detailed information about the proposal steps, including:
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
#'   \item{\code{unif_scaler}}{A numeric scalar, the inverse probability of sampling from the steps part of the proposal (\eqn{\frac{1}{p(lower < x < upper)}}). Used for scaling uniform random values.}
#'   \item{\code{lt_properties}}{A numeric vector of 5 values required for Adaptive Rejection Sampling (ARS) in the left tail.}
#'   \item{\code{rt_properties}}{A numeric vector of 6 values required for ARS in the right tail.}
#'   \item{\code{alpha}}{A numeric scalar representing the uniform step area.}
#'   \item{\code{tails_method}}{A string, either \code{"ARS"} (Adaptive Rejection Sampling) or \code{"IT"} (Inverse Transform), indicating the sampling method for the tails.}
#'   \item{\code{proposal_bounds}}{A numeric vector specifying the left and right bounds of the target density.}
#'   \item{\code{cnum}}{An integer representing the cache number of the created proposal in memory.}
#'   \item{\code{symmetric}}{A numeric scalar indicating the symmetry point of the proposal, or \code{NULL} if not symmetric.}
#'   \item{\code{f_params}}{A list of parameters for the target density that the proposal is designed for.}
#'   \describe{
#'   \item{\code{mu}}{is the location parameter (location of the distribution).}
#'   \item{\code{b}}{is the scale parameter, which controls the spread of the distribution (\code{b > 0}).}
#'   }
#'   \item{\code{is_symmetric}}{A logical value indicating whether the proposal is symmetric.}
#'   \item{\code{proposal_type}}{A string indicating the type of the generated proposal:
#'   \describe{
#'     \item{\code{"scaled"}}{The proposal is "scalable" and standardized with \code{mu = 0} and \code{b = 1}. This is used when parameters \code{mu} and \code{b} are either \code{NULL} or not provided. Scalable proposals are compatible with \code{\link{srlaplace}}.}
#'     \item{\code{"custom"}}{The proposal is "custom" when either \code{mu} or \code{b} is provided. Custom proposals are compatible with \code{\link{srlaplace_custom}}.}
#'   }}
#'   \item{\code{target_function_area}}{A numeric scalar estimating the area of the target distribution.}
#'   \item{\code{dens_func}}{A string containing the hardcoded density function.}
#'   \item{\code{density_name}}{A string specifying the name of the target density distribution.}
#'   \item{\code{lock}}{An identifier used for saving and loading the proposal from disk.}
#' }
#'
#' @seealso
#' \code{\link{srlaplace}}: Function to sample from a scalable proposal generated by \code{srlaplace_optimize()}.
#' \code{\link{srlaplace_custom}}: Function to sample from a custom proposal tailored to user specifications.
#'
#'
#' @examples
#' # Generate scalable proposal that with mu = 0 and b = 1, that has 4096 steps
#' scalable_proposal <- srlaplace_optimize(steps = 4096)
#'
#' # Generate custom proposal that with mu = 2 and b = 1
#' scalable_proposal <- srlaplace_optimize(mu = 2, b = 1)
#'
#' @export
srlaplace_optimize <- function(
    mu = NULL,
    b = NULL,
    xl = NULL,
    xr = NULL,
    steps = 4091,
    proposal_range = NULL,
    theta = NULL,
    target_sample_size = 1000,
    verbose = FALSE,
    symmetric = FALSE
) {

  dist_name <- "srlaplace"

  dendata <- built_in_proposals[[dist_name]]

  f_params <- list(mu = mu, b = b)

  if (dendata$scalable) {

    isnull <- sapply(f_params, is.null)

    if (all(isnull)) {
      cnum <- dendata$c_num
      proposal_type <- "scaled"
    } else {
      cnum <- dendata$c_num + 1
      proposal_type <- "custom"
    }

    f_params <- ifelse(isnull, dendata$std_params, f_params)

  } else {

    cnum <- dendata$c_num + 1
    proposal_type <- "custom"
  }

  modes <- dendata$set_modes(f_params$mu)

  if (symmetric) {
    symmetric <- modes
  } else {
    symmetric <- NULL
  }

  f <- dendata$create_f(f_params$mu, f_params$b)

  check_proposal_opt_criteria(symmetric, cnum, dendata)

  proposal_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 proposal_range, theta, target_sample_size,
                 proposal_type, symmetric, cnum, verbose)

}

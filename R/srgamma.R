#' Sampling from Gamma Distribution
#' @rdname srgamma_custom
#'
#' @description
#' The \code{srgamma_custom()} function generates random samples from a Gamma distribution using the STORS algorithm.
#' It employs an optimized proposal distribution around the mode and Adaptive Rejection Sampling (ARS) for the tails.
#'
#' @details
#'
#' The Gamma Distribution
#' @details
#' The Gamma distribution has the probability density function (PDF):
#' \deqn{f(x | \alpha, \beta) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha - 1} \exp(-\beta x), \quad x \geq 0,}
#' where:
#' \describe{
#'   \item{\eqn{\alpha}}{is the shape parameter (\eqn{\alpha > 0}), which determines the shape of the distribution.}
#'   \item{\eqn{\beta}}{is the rate parameter (\eqn{\beta > 0}), which determines the rate of decay.}
#' }
#' The Gamma distribution is widely used in statistics, particularly in Bayesian inference and modeling waiting times.
#'
#' This function samples from a proposal constructed using \code{\link{srgamma_optimize}}, employing the STORS algorithm.
#'
#' By default, \code{srgamma_custom()} samples from the standard Gamma distribution with \code{shape = 1} and \code{rate = 1}.
#' The proposal distribution is pre-optimized at package load time using \code{srgamma_optimize()} with
#' \code{steps = 4091}, creating a scalable proposal centered around the mode.
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is overwritten in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing random samples from the Gamma distribution.
#' The \code{shape} and \code{rate} parameters are specified during the optimization process using \code{srgamma_optimize()}.
#'
#' \bold{NOTE:} When the \code{x} parameter is specified, it is updated in-place with the simulation for performance reasons.
#'
#' @note
#' This function is not scalable. Therefore, only the \code{srgamma_custom()} version is available, which requires the proposal to be pre-optimized using \code{srgamma_optimize()} before calling this function.
#'
#' @seealso
#' \code{\link{srgamma_optimize}} to optimize the custom proposal.
#'
#' @examples
#' # Generate 10 samples from Gamma Distribution
#' samples <- srgamma_custom(10)
#' print(samples)
#'
#' # Generate 10 samples using a pre-allocated vector
#' x <- numeric(10)
#' srgamma_custom(10, x = x)
#' print(x)
#'
#' @export
srgamma_custom <- function(n = 1, x = NULL) {
  .Call(C_srgamma_custom_check, n, x)
}


#' Optimizing Gamma Distribution proposal
#' @description
#' The \code{srgamma_optimize()} function generates an optimized proposal for a targeted Gamma distribution.
#' The proposal can be customized and adjusted based on various options provided by the user.
#'
#' @details
#' When \code{srgamma_optimize()} is explicitly called:
#' \itemize{
#'   \item A proposal is created and cached. If no parameters are provided, a standard proposal is created with \code{shape = 1} and \code{rate = 1}.
#'   \item Providing \code{shape} and \code{rate} creates a custom proposal, which is cached for use with \code{srgamma_custom()}.
#'   \item The optimization process can be controlled via parameters such as \code{steps}, \code{proposal_range}, or
#'   \code{theta}. If no parameters are provided, the proposal is optimized via brute force based on the
#'   \code{target_sample_size}.
#' }
#'
#' @param shape (optional) Numeric. The shape parameter (\eqn{\alpha}) of the Gamma distribution. Defaults to \code{1}.
#' @param rate (optional) Numeric. The rate parameter (\eqn{\beta}) of the Gamma distribution. Defaults to \code{1}.
#' @param scale (optional) Numeric. The scale parameter of the Gamma distribution. Defaults to \code{1}.
#' @param xl Numeric. Left truncation bound for the target distribution. Defaults to \code{0}, representing no left truncation.
#' @param xr Numeric. Right truncation bound for the target distribution. Defaults to \code{Inf}, representing no right truncation.
#' @param steps (optional) Integer. Desired number of steps in the proposal. Defaults to \code{4091}.
#' @param proposal_range (optional) Numeric vector. Specifies the range for optimizing the steps part of the proposal. Defaults to \code{NULL}, indicating automatic range selection.
#' @param theta (optional) Numeric. A parameter for proposal optimization. Defaults to \code{NULL}.
#' @param target_sample_size (optional) Integer. Target sample size for proposal optimization. Defaults to \code{1000}.
#' @param verbose Boolean. If \code{TRUE}, detailed optimization information, including areas and steps, will be displayed. Defaults to \code{FALSE}.
#'
#' @return
#' A list containing the optimized proposal and related parameters for the specified Gamma distribution. The proposal is also cached for internal use.
#' \describe{
#'   \item{\code{data}}{Detailed information about the proposal steps, including \code{x}, \code{s_upper}, \code{p_a}, and \code{s_upper_lower}.}
#'   \item{\code{areas}}{The areas under the left tail, steps, and right tail of the proposal distribution.}
#'   \item{\code{steps_number}}{The number of steps in the proposal.}
#'   \item{\code{f_params}}{The parameters (\code{shape} and \code{rate}) of the Gamma distribution.}
#' }
#'
#' @seealso
#' \code{\link{srgamma_custom}}: Function to sample from a custom proposal generated by \code{srgamma_optimize()}.
#'
#' @examples
#' # Generate a standard proposal with shape = 1 and rate = 1
#' standard_proposal <- srgamma_optimize()
#'
#' # Generate a custom proposal with shape = 2 and rate = 3
#' custom_proposal <- srgamma_optimize(shape = 2, rate = 3)
#'
#' @export
srgamma_optimize <- function(shape = NULL,
                            rate = NULL,
                            scale = NULL,
                            xl = NULL,
                            xr = NULL,
                            steps = 4091,
                            proposal_range = NULL,
                            theta = NULL,
                            target_sample_size = 1000,
                            verbose = FALSE) {

  dist_name <- "srgamma"

  symmetric <- NULL

  dendata <- built_in_proposals[[dist_name]]

  cnum <- dendata$c_num + 1

  proposal_type <- "custom"

  if (!is.null(rate)) scale <- 1 / rate

  f_params <- list(shape = shape, scale = scale) # F L

  f_params <- ifelse(sapply(f_params, is.null), dendata$std_params, f_params)

  modes <- dendata$set_modes(f_params$shape, f_params$scale)

  f <- dendata$create_f(shape = f_params$shape, scale = f_params$scale)

  proposal_optimizer(
    dendata,
    dist_name,
    xl,
    xr,
    f,
    modes,
    f_params,
    steps,
    proposal_range,
    theta,
    target_sample_size,
    proposal_type,
    symmetric,
    cnum,
    verbose
  )

}

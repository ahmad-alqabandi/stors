#' Sampling from Pareto Distribution
#' @rdname srpareto_custom
#'
#' @description
#' The \code{srpareto_custom()} function generates random samples from a Pareto distribution using the STORS algorithm.
#' It employs an optimized proposal distribution around the mode and Inverse Transform (IT) method for the tails.
#'
#' @details
#'
#' The Pareto Distribution
#'
#' @details
#' The Pareto distribution has the probability density function (PDF):
#' \deqn{f(x | \alpha, \sigma) = \frac{\alpha \sigma^\alpha}{x^{\alpha + 1}}, \quad x \geq \sigma,}
#' where:
#' \describe{
#'   \item{\eqn{\alpha}}{is the shape parameter (\eqn{\alpha > 0}), which determines the tail heaviness of the distribution.}
#'   \item{\eqn{\sigma}}{is the scale parameter (\eqn{\sigma > 0}), which determines the minimum possible value of \eqn{x}.}
#' }
#' The Pareto distribution is widely used in modelling phenomena with heavy tails, such as wealth distribution, insurance losses, and natural events.
#'
#' This function samples from a proposal constructed using \code{\link{srpareto_optimize}}, employing the STORS algorithm.
#'
#' By default, \code{srpareto_custom()} samples from the standard Pareto distribution with \code{shape = 1} and \code{rate = 1}.
#' The proposal distribution is pre-optimized at package load time using \code{srpareto_optimize()} with
#' \code{steps = 4091}, creating a scalable proposal centred around the mode.
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is overwritten in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing random samples from the Pareto distribution.
#' The \code{shape} and \code{scale} parameters are specified during the optimization process using \code{srpareto_optimize()}.
#'
#' \bold{NOTE:} When the \code{x} parameter is specified, it is updated in-place with the simulation for performance reasons.
#'
#' @note
#' This function is not scalable. Therefore, only the \code{srpareto_custom()} version is available, which requires the proposal to be pre-optimized using \code{srpareto_optimize()} before calling this function.
#'
#' @seealso
#' \code{\link{srpareto_optimize}} to optimize the custom proposal.
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



#' Optimizing Pareto Distribution proposal
#' @description
#' The \code{srpareto_optimize()} function generates an optimized proposal for a targeted Pareto Distribution.
#'  The proposal can be customized and adjusted based on various options provided by the user.
#'
#'
#' @details
#'When \code{srpareto_optimize()} is explicitly called:
#'\itemize{
#'  \item A proposal is created and cached. If no parameters are provided, a standard proposal is created with \code{rate = 1}.
#'  \item Providing \code{rate} creates a custom proposal, which is cached for use with \code{srpareto_custom()}.
#'  \item The optimization process can be controlled via parameters such as \code{steps}, \code{proposal_range}, or
#'   \code{theta}. If no parameters are provided, the proposal is optimized via brute force based on the.
#'   \code{target_sample_size}.
#'}
#'
#' @param shape (optional) Numeric. shape parameter of the Pareto Distribution. Defaults to \code{NULL}, which implies a scalable proposal with \code{shape = 1}.
#' @param scale (optional) Numeric. scale parameter of the Pareto Distribution. Defaults to \code{NULL}, which implies a scalable proposal with \code{scale = 1}.
#' @param xl Numeric. Left truncation bound for the target distribution. Defaults to \code{-Inf}, representing no left truncation.
#' @param xr Numeric. Right truncation bound for the target distribution. Defaults to \code{Inf}, representing no right truncation.
#' @param steps (optional) Integer. Desired number of steps in the proposal. Defaults to \code{NULL}, which means the number of steps is determined automatically during optimization.
#' @param proposal_range (optional) Numeric vector. Specifies the range for optimizing the steps part of the proposal. Defaults to \code{NULL}, indicating automatic range selection.
#' @param theta Numeric. A parameter for proposal optimization. Defaults to 0.1.
#' @param target_sample_size (optional) Integer. Target sample size for proposal optimization. Defaults to \code{1000}.
#' @param verbose Boolean. If \code{TRUE}, detailed optimization information, including areas and steps, will be displayed. Defaults to \code{FALSE}.
#'
#'
#' @return
#' The user does not need to store the returned value, because the package internally cashes the proposal. However, we explain here the full returned proposal for advanced users.
#'
#' A list containing the optimized proposal and related parameters for the specified built-in distribution:
#'#' \describe{
#'   \item{\code{data}}{Detailed information about the proposal steps, including \code{x}, \code{s_upper}, \code{p_a}, and \code{s_upper_lower}.}
#'   \item{\code{areas}}{The areas under the left tail, steps, and right tail of the proposal distribution.}
#'   \item{\code{steps_number}}{The number of steps in the proposal.}
#'   \item{\code{f_params}}{The parameters (\code{scale} and \code{shape}) of the Beta distribution.}
#' }
#'
#' @seealso
#' \code{\link{srpareto_custom}}: Function to sample from a custom proposal tailored to user specifications.
#'
#'
#' @examples
#' # Generate scalable proposal that with rate = 1, that has 4096 steps
#' scalable_proposal <- srpareto_optimize(steps = 4096)
#'
#' # Generate custom proposal that with scale = 4
#' scalable_proposal <- srpareto_optimize(scale = 4)
#'
#' @export
srpareto_optimize <- function(scale = NULL,
                              shape = NULL,
                              xl = NULL,
                              xr = NULL,
                              steps = 4091,
                              proposal_range = NULL,
                              theta = 0.1,
                              target_sample_size = 1000,
                              verbose = FALSE) {
  dist_name <- "srpareto"

  dendata <- built_in_proposals[[dist_name]]

  f_params <- list(scale = scale, shape = shape)

  cnum <- dendata$c_num + 1

  proposal_type <- "custom"

  f_params <- ifelse(sapply(f_params, is.null), dendata$std_params, f_params)

  modes <- dendata$set_modes(f_params$scale)

  symmetric <- NULL

  f <- dendata$create_f(f_params$scale, f_params$shape)

  dendata$lower <- f_params$scale

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

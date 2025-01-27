#' Sampling from Chi-squared Distribution
#' @rdname srchisq_custom
#'
#' @description
#' The \code{srchisq_custom()} function generates random samples from a Chi-squared Distribution using the STORS algorithm.
#' It employs an optimized proposal distribution around the mode and Adaptive Rejection Sampling (ARS) for the tails.
#'
#' @details
#'
#' The Chi-squared Distribution
#' @details
#' The Chi-squared distribution has the probability density function (PDF):
#' \deqn{f(x | k) = \frac{1}{2^{k/2} \Gamma(k/2)} x^{(k/2) - 1} \exp(-x/2), \quad x \geq 0,}
#' where:
#' \describe{
#'   \item{\eqn{k}}{ is the degrees of freedom (\eqn{k > 0}), which determines the shape of the distribution.}
#' }
#' The Chi-squared distribution is widely used in hypothesis testing and constructing confidence intervals, particularly in the context of variance estimation.
#'
#' this function is sampling from proposal that has been constructed using \code{\link{srchisq_optimize}}, using the STORS algorithm.
#'
#' By default, \code{srchisq_custom()} samples from Chi-squared Distribution \code{df = 2}.
#' The proposal distribution is pre-optimized at package load time using \code{srchisq_optimize()} with
#' \code{steps = 4091}, creating a scalable proposal centered around the mode.
#'
#'
#' @param n Integer, length 1. Number of samples to draw.
#' @param x (optional) Numeric vector of length \eqn{n}. If provided, this vector is over
#' written in place to avoid any memory allocation.
#'
#' @return
#' A numeric vector of length \code{n} containing random samples from the Chi-squared distribution.
#' The degrees of freedom (\code{df}) for the distribution are specified during the optimization process using \code{srchisq_optimize()}.
#' \bold{NOTE:} When the \code{x} parameter is specified, it is updated in-place with the simulation for performance reasons.
#'
#' @note
#' This function is not scalable. Therefore, only the \code{srchisq_custom()} version is available, which requires the proposal to be pre-optimized using \code{srchisq_optimize()} before calling this function.
#'
#' @seealso
#' \code{\link{srchisq_optimize}} to optimize the custom proposal.
#'
#' @examples
#'
#' # Genedf 10 samples from Chi-squared Distribution
#' samples <- srchisq_custom(10)
#' print(samples)
#'
#' # Genedf 10 samples using a pre-allocated vector
#' x <- numeric(10)
#' srchisq_custom(10, x = x)
#' print(x)
#'
#'
#' @export
srchisq_custom <- function(n = 1, x = NULL) {
  .Call(C_srchisq_custom_check, n, x)
}




#' Optimizing Chi-squared Distribution proposal
#' @description
#' The \code{srchisq_optimize()} function generates an optimized proposal for a targeted Chi-squared Distribution.
#'  The proposal can be customized and adjusted based on various options provided by the user.
#'
#'
#' @details
#'When \code{srchisq_optimize()} is explicitly called:
#' \itemize{
#'  \item A proposal is created and cached. If no parameters are provided, a standard proposal is created with \code{df = 2}.
#'  \item Providing \code{df} creates a custom proposal, which is cached for use with \code{srchisq_custom()}.
#'  \item The optimization process can be controlled via parameters such as \code{steps}, \code{proposal_range}, or
#'   \code{theta}. If no parameters are provided, the proposal is optimized via brute force based on the.
#'   \code{target_sample_size}.
#'}
#'
#' @param df (optional) Numeric. degrees of freedom parameter of the Chi-squared Distribution. Defaults to \code{NULL}, which implies proposal with \code{df = 2}.
#' @param xl Numeric. Left truncation bound for the target distribution. Defaults to \code{-Inf}, representing no left truncation.
#' @param xr Numeric. Right truncation bound for the target distribution. Defaults to \code{Inf}, representing no right truncation.
#' @param steps (optional) Integer. Desired number of steps in the proposal. Defaults to \code{NULL}, which means the number of steps is determined automatically during optimization.
#' @param proposal_range (optional) Numeric vector. Specifies the range for optimizing the steps part of the proposal. Defaults to \code{NULL}, indicating automatic range selection.
#' @param theta (optional) Numeric. A parameter for proposal optimization. Defaults to \code{NULL}.
#' @param target_sample_size (optional) Integer. Target sample size for proposal optimization. Defaults to \code{1000}.
#' @param verbose Boolean. If \code{TRUE}, detailed optimization information, including areas and steps, will be displayed. Defaults to \code{FALSE}.
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
#'     \item{\code{df}}{ the df of the target distribution.}
#'   }
#'   \item{\code{is_symmetric}}{A logical value indicating whether the proposal is symmetric.}
#'   \item{\code{proposal_type}}{A string indicating the type of the genedfd proposal:
#'   \describe{
#'     \item{\code{"custom"}}{The proposal is "custom" when \code{df} is provided. Custom proposals are compatible with \code{\link{srchisq_custom}}.}
#'   }}
#'   \item{\code{target_function_area}}{A numeric scalar estimating the area of the target distribution.}
#'   \item{\code{dens_func}}{A string containing the hardcoded density function.}
#'   \item{\code{density_name}}{A string specifying the name of the target density distribution.}
#'   \item{\code{lock}}{An identifier used for saving and loading the proposal from disk.}
#' }
#'
#' @seealso
#' \code{\link{srchisq_custom}}: Function to sample from a custom proposal genedfd by \code{srchisq_optimize()}.
#'
#'
#' @examples
#' # Genedf custom proposal that with df = 2, that has 4096 steps
#' scalable_proposal <- srchisq_optimize(steps = 4096)
#'
#' # Genedf custom proposal that with df = 4
#' scalable_proposal <- srchisq_optimize(df = 4)
#'
#' @export
srchisq_optimize <- function(df = 2,
                            xl = NULL,
                            xr = NULL,
                            steps = 4091,
                            proposal_range = NULL,
                            theta = NULL,
                            target_sample_size = 1000,
                            verbose = FALSE) {

  dist_name <- "srchisq"

  symmetric <- NULL

  dendata <- built_in_proposals[[dist_name]]

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

  if (df == 1) {
    message(
      "proposal building is not available for df = 1. You can square the result of srnorm() to sample from this distribution."
    )
    return()
  }

  f_params <- list(df = df)

  modes <- dendata$set_modes(df)

  f <- dendata$create_f(df)

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

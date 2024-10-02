#' Sampling from Exponential Distribution
#' @rdname srexp
#' @order 1
#'
#' @description
#' Sampling from the exponential distribution using Stors.
#'
#' @details
#' The function \code{srexp()} is used for sampling from a standard exponential distribution with a rate parameter of 1.
#'
#' Exponential distribution has the density:
#' \deqn{ f(x; \lambda) = \lambda e^{-\lambda x} }
#' Where \eqn{\lambda} is the rate parameter.
#'
#' @param n Integer sample size.
#'
#' @return
#' \code{srexp()} returns a sample of size \code{n} from a standard exponential distribution.
#'
#' @examples
#' 
#' # Optimize the grid for Exponential distribution
#' grid_optimizer("srexp")
#' 
#' # Generating Samples from a Standard Exponential Distribution
#' # This example illustrates how to generate 10 samples,
#' # from a standard exponential distribution.
#' samples <- srexp(10)
#' print(samples)
#'
#' @export
srexp <- function(n) {
  .Call(C_srexp, n)
}


#' @rdname srexp
#' @order 2
#'
#' @description
#' Sampling from any exponential distribution by specifying the scale parameter.
#'
#' @param n Integer sample size.
#' @param scale Scalar scale parameter.
#'
#' @return
#' \code{srexp_scaled()} returns a sample of size \code{n} from an exponential distribution with the specified scale parameter.
#'
#' @examples
#' # Generating Samples from an Exponential Distribution with a Specific Scale
#' # This example demonstrates how to generate 10 samples,
#' # from an exponential distribution with a scale parameter of 2.
#' samples <- srexp_scaled(n = 10, scale = 2)
#' print(samples)
#'
#' @export
srexp_scaled <- function(n, scale = 1) {
  .Call(C_srexp, n) * scale
}



#' @rdname srexp
#' @order 3
#'
#' @description
#' Creating a sampling function for a truncated exponential distribution.
#'
#' @details
#' \code{srexp_truncate()} is used for sampling from a standard exponential distribution truncated within specified bounds.
#' The function validates the truncation bounds before creating the sampling function.
#'
#' @param xl Lower bound for truncation.
#' @param xr Upper bound for truncation.
#'
#' @return
#' \code{srexp_truncate()} returns a function that generates \code{n} samples from a truncated exponential distribution between \code{xl} and \code{xr}.
#'
#' @examples
#' # Generating Samples from a Truncated Exponential Distribution
#' # This example demonstrates how to generate 100 samples,
#' # from an exponential distribution truncated in the range [1, 3].
#' exp_trunc <- srexp_truncate(xl = 1, xr = 3)
#' samples <- exp_trunc(100)
#' hist(samples, main = "Histogram of Truncated Exponential Samples", xlab = "Value", breaks = 20)
#'
#' @export
srexp_truncate = function(xl, xr){
  
  truncate_error_checking(xl, xr, pbgrids$srexp$lb, pbgrids$srexp$rb)
  
  Upper_cumsum = .Call(C_srexp_trunc_nav, xl, xr)
  
  function_string <- paste0("function(n) { .Call(C_srexp_trunc, n, ", paste0(xl), ", ", paste0(xr), ", ", paste0(Upper_cumsum[1]), ", ", paste0(Upper_cumsum[2]), ", ", paste0(as.integer(Upper_cumsum[3])), ", ", paste0(as.integer(Upper_cumsum[4])), ")}")
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  
}

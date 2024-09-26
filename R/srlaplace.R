

#' Sampling from Laplace Distribution
#' @rdname srlaplace
#' @order 1
#'
#' @description
#' Sampling from the Laplace distribution using Stors.
#'
#' @details
#' The function `srlaplace()` is used for sampling from a standard Laplace distribution (location = 0, scale = 1).
#'
#' Laplace distribution has the density:
#' \deqn{ f(x | \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right) }
#' Where \eqn{\mu} is the location parameter and \eqn{b} is the scale parameter.
#'
#' @param n Integer sample size.
#'
#' @return
#' `srlaplace()` returns a sample of size `n` from a standard Laplace distribution.
#'
#' @examples
#' 
#' # Optimize the grid for the Laplace distribution
#' grid_optimizer("laplace")
#' 
#' # Generating Samples from a Standard Laplace Distribution
#' # This example illustrates how to generate 10 samples from a standard Laplace distribution.
#' samples <- laplace(10)
#' print(samples)
#' 
#'
#' @export
srlaplace <- function(n) {
  .Call(C_laplace, n)
}



#' @rdname srlaplace
#' @order 2
#'
#' @description
#' Sampling from any Laplace distribution specifying the location and scale.
#'
#' @param mu Scalar location parameter.
#' @param b Scalar scale parameter.
#'
#' @return
#' `srlaplace_scaled()` returns a sample of size `n` from a Laplace distribution with location \eqn{\mu} and scale \eqn{b}.
#'
#' @examples
#' # Generating Samples from a Laplace Distribution with Specific Location and Scale
#' # This example demonstrates how to generate 10 samples from a Laplace distribution with location = 2 and scale = 3.
#' samples <- srlaplace_scaled(n = 10, mu = 2, b = 3)
#' print(samples)
#'
#' @export
srlaplace_scaled <- function(n, mu = 0, b = 1) {
  .Call(C_laplace, n) * b + mu
}



#' @rdname srlaplace
#' @order 3
#'
#' @description
#' Creating a sampling function for a truncated Laplace distribution.
#'
#' @details
#' `srlaplace_truncate()` is used for sampling from a standard Laplace distribution truncated within specified bounds. 
#' It is beneficial when the area of interest in a Laplace distribution is limited to a specific range. The function validates 
#' the truncation bounds before creating the sampling function.
#'
#' @param xl Lower bound for truncation.
#' @param xr Upper bound for truncation.
#'
#' @return
#' `srlaplace_truncate()` returns a function that generates `n` samples from a truncated Laplace distribution between `xl` and `xr`.
#'
#' @examples
#' # Generating Samples from a Truncated Laplace Distribution
#' # This example demonstrates how to generate 100 samples from a Laplace distribution truncated in the range [1, 3].
#' laplace_trunc <- srlaplace_truncate(xl = 1, xr = 3)
#' samples <- laplace_trunc(100)
#' hist(samples, main = "Histogram of Truncated Laplace Samples", xlab = "Value", breaks = 20)
#'
#' @export
srlaplace_truncate = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >=  pbgrids$laplace$lb,
    "xr must be smaller than the density upper bound" = xr <= pbgrids$laplace$rb
  )
  
  Upper_cumsum = .Call(C_laplace_trunc_nav, xl, xr)
  
  function_string <- paste0("function(n) { .Call(C_laplace_trunc, n, ", paste0(xl), ", ", paste0(xr), ", ", paste0(Upper_cumsum[1]), ", ", paste0(Upper_cumsum[2]), ", ", paste0(as.integer(Upper_cumsum[3])), ", ", paste0(as.integer(Upper_cumsum[4])), ")}")
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
}



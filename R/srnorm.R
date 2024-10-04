#' Sampling from Normal Distribution
#' @rdname srnorm
#' @order 1
#'
#' @description
#' Sampling from the Normal distribution using Stors.
#'
#' @details
#' The function \code{srnorm()} is used for sampling from a standard Normal distribution (mean = 0, standard deviation = 1).
#' 
#' 
#'
#' Normal distribution has the density:
#' 
#' \deqn{ f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mu)^2}{2\sigma^2}} }
#' 
#' Where \eqn{\mu} is the mean and \eqn{\sigma} is the standard deviation.
#'
#'
#' @param n Integer sample size.
#'
#' @return
#' \code{srnorm()} returns a sample of size \code{n} from a standard normal distribution.
#'
#' @examples
#' # Generating Samples from a Standard Normal Distribution
#' # This example illustrates how to generate 10 samples from a standard normal distribution.
#' # It first optimizes the grid for sampling using \code{grid_optimizer},
#' # and then generates samples using \code{srnorm}.
#'
#' # Optimize the grid for the standard normal distribution
#' grid_optimizer("srnorm")
#'
#' # Generate and print 10 samples from the standard normal distribution
#' samples <- srnorm(10)
#' print(samples)
#'
#'
#' @export
srnorm <- function(n) {
  .Call(C_srnorm, n)
}

#' @rdname srnorm
#' @order 2
#'
#' @param mu Scalar mean.
#' @param sd Scalar standard deviation.
#' 
#' 
#' @details
#' \code{srnorm_scaled()} allows sampling from any Normal distribution by specifying the mean and standard deviation.
#' The separation of these functions enhances performance, as the Stors algorithm is highly efficient, and even simple arithmetic can impact its speed.
#'
#' @return
#' \code{srnorm_scaled()} returns a sample of size \code{n} from a normal distribution with mean \eqn{\mu} and standard deviation \eqn{\sigma}.
#'
#' @examples
#' # Generating Samples from a Normal Distribution with Specific Mean and Standard Deviation
#' # This example demonstrates how to generate 10,
#' # samples from a normal distribution with a mean of 4 and a standard deviation of 2.
#'
#' samples <- srnorm_scaled(n = 10, mu = 4, sd = 2)
#' print(samples)
#' 
#' @export
srnorm_scaled <- function(n, mu = 0, sd = 1) {
  .Call(C_srnorm, n) * sd + mu
}

#' @rdname srnorm
#' @order 3
#'
#' @param xl Lower bound for truncation.
#' @param xr Upper bound for truncation.
#' 
#' 
#' @details
#' \code{srnorm_truncate()}, this function allows sampling from a standard normal distribution that is truncated within specified bounds.
#'  It is particularly useful when the area of interest in a normal distribution is limited to a specific range.
#'  The function first validates the truncation bounds to ensure they are within the allowable range of the distribution and then creates a tailored sampling function based on these bounds.
#' 
#' @return
#' \code{srnorm_truncate()} returns a function that, when called with a sample size \code{n}, generates \code{n} samples from a normal distribution truncated between \code{xl} and \code{xr}.
#'
#' @examples
#' # Generating Samples from a Truncated Standard Normal Distribution
#' # This example demonstrates how to generate 100,
#' # samples from a standard normal distribution truncated in the range [-2, 2].
#'
#' # Create the truncated sampling function
#' norm_trunc <- srnorm_truncate(xl = -2, xr = 2)
#'
#' # Generate 100 samples
#' sample <- norm_trunc(100)
#'
#' # Plot a histogram of the samples
#' hist(sample, main = "Histogram of Truncated Normal Samples", xlab = "Value", breaks = 20)
#' 
#' @export
srnorm_truncate <- function(xl = -Inf, xr = Inf) {
  
  
  res <- truncate_error_checking(xl, xr, pbgrids$srnorm)
  xl <- res[[1]]
  xr <- res[[2]]
  
  Upper_cumsum = .Call(C_srnorm_trunc_nav, xl, xr)
  
  function_string <- paste0("function(n) { .Call(C_srnorm_trunc, n, ", paste0(xl), ", ", paste0(xr), ", ", paste0(Upper_cumsum[1]), ", ", paste0(Upper_cumsum[2]), ", ", paste0(as.integer(Upper_cumsum[3])), ", ", paste0(as.integer(Upper_cumsum[4])), ")}")
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
}

#' Sampling from Normal Distribution
#' @rdname srgamma
#' @order 1
#'
#' @description
#' Sampling from the Normal distribution using Stors.
#'
#' @details
#' The function \code{srgamma()} is used for sampling from a standard Normal distribution (mean = 0, standard deviation = 1).
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
#' \code{srgamma()} returns a sample of size \code{n} from a standard normal distribution.
#'
#' @examples
#' # Generating Samples from a Standard Normal Distribution
#' # This example illustrates how to generate 10 samples from a standard normal distribution.
#' # It first optimizes the grid for sampling using \code{grid_optimizer},
#' # and then generates samples using \code{srgamma}.
#'
#' # Optimize the grid for the standard normal distribution
#' grid_optimizer("srgamma")
#'
#' # Generate and print 10 samples from the standard normal distribution
#' samples <- srgamma(10)
#' print(samples)
#'
#'
#' @export
srgamma <- function(n) {
  .Call(C_srgamma, n)
}


#' @rdname srgamma
#' @order 2
#'
#' @param xl Lower bound for truncation.
#' @param xr Upper bound for truncation.
#' 
#' 
#' @details
#' \code{srgamma_truncate()}, this function allows sampling from a standard normal distribution that is truncated within specified bounds.
#'  It is particularly useful when the area of interest in a normal distribution is limited to a specific range.
#'  The function first validates the truncation bounds to ensure they are within the allowable range of the distribution and then creates a tailored sampling function based on these bounds.
#' 
#' @return
#' \code{srgamma_truncate()} returns a function that, when called with a sample size \code{n}, generates \code{n} samples from a normal distribution truncated between \code{xl} and \code{xr}.
#'
#' @examples
#' # Generating Samples from a Truncated Standard Normal Distribution
#' # This example demonstrates how to generate 100,
#' # samples from a standard normal distribution truncated in the range [-2, 2].
#'
#' # Create the truncated sampling function
#' norm_trunc <- srgamma_truncate(xl = -2, xr = 2)
#'
#' # Generate 100 samples
#' sample <- norm_trunc(100)
#'
#' # Plot a histogram of the samples
#' hist(sample, main = "Histogram of Truncated Normal Samples", xlab = "Value", breaks = 20)
#' 
#' @export
srgamma_truncate <- function(xl = 0, xr = Inf) {
  
  
  res <- truncate_error_checking(xl, xr, pbgrids$srgamma)
  xl <- res$xl; xr <- res$xr
  
  Upper_cumsum = .Call(C_srgamma_trunc_nav, xl, xr)
  
  stopifnot(
    "xl is has a CDF close to 1" = (Upper_cumsum[1] != 1),
    "xr is has a CDF close to 0" = (Upper_cumsum[2] != 0)
  )
  
  function_string <- paste0("function(n) { .Call(C_srgamma_trunc, n, ", paste0(xl), ", ", paste0(xr), ", ", paste0(Upper_cumsum[1]),
                            ", ", paste0(Upper_cumsum[2]), ", ", paste0(as.integer(Upper_cumsum[3])), ", ", paste0(as.integer(Upper_cumsum[4])),")}")
  
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
}



#' @export
srgamma_optimize = function(
  shape = 1,
  rate = 1,
  scale = 1/rate,
  steps = 4091,
  grid_range = NULL,
  theta = NULL,
  target_sample_size = 1000,
  verbose = FALSE
) {
  
  density_name <- 'srgamma'
  
  dendata <- pbgrids[[density_name]]
  
  f_params <- c(shape, scale) # F L
  
  modes <- dendata$set_modes(shape, scale)

  f <- dendata$create_f(shape = shape, scale = scale)
  
  if( identical(dendata$tails_method,"ARS") ){
    
    h <- function(x)
      log(f(x))
    
    h_prime <- stors_prime(modes, h)
    
  }else{
    cdf <- dendata$create_cdf(mu, sd)
  }  
  
  grid_optimizer(dendata, density_name, f, cdf, h,
                 h_prime, modes, f_params, steps,
                 grid_range, theta, target_sample_size, verbose)
  
}



#' Sampling from Laplace Distribution
#' @rdname srlaplace
#' @order 1
#'
#' @description
#' Sampling from the Laplace distribution using Stors.
#'
#' @details
#' The function \code{srlaplace()} is used for sampling from a standard Laplace distribution (location = 0, scale = 1).
#'
#' Laplace distribution has the density:
#' \deqn{ f(x | \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right) }
#' Where \eqn{\mu} is the location parameter and \eqn{b} is the scale parameter.
#'
#' @param n Integer sample size.
#'
#' @return
#' \code{srlaplace()} returns a sample of size \code{n} from a standard Laplace distribution.
#'
#' @examples
#' 
#' # Optimize the grid for the Laplace distribution
#' grid_optimizer("srlaplace")
#' 
#' # Generating Samples from a Standard Laplace Distribution
#' # This example illustrates how to generate 10 samples from a standard Laplace distribution.
#' samples <- srlaplace(10)
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
#' \code{srlaplace_scaled()} returns a sample of size \code{n} from a Laplace distribution with location \eqn{\mu} and scale \eqn{b}.
#'
#' @examples
#' # This example demonstrates how to generate 10,
#' # Generating Samples from a Laplace Distribution with Specific Location and Scale
#' # samples from a Laplace distribution with location = 2 and scale = 3.
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
#' \code{srlaplace_truncate()} is used for sampling from a standard Laplace distribution truncated within specified bounds. 
#' It is beneficial when the area of interest in a Laplace distribution is limited to a specific range. The function validates 
#' the truncation bounds before creating the sampling function.
#'
#' @param xl Lower bound for truncation.
#' @param xr Upper bound for truncation.
#'
#' @return
#' \code{srlaplace_truncate()} returns a function that generates \code{n} samples from a truncated Laplace distribution between \code{xl} and \code{xr}.
#'
#' @examples
#' # Generating Samples from a Truncated Laplace Distribution
#' # This example demonstrates how to generate 100 samples,
#' # from a Laplace distribution truncated in the range [1, 3].
#' laplace_trunc <- srlaplace_truncate(xl = 1, xr = 3)
#' samples <- laplace_trunc(100)
#' hist(samples, main = "Histogram of Truncated Laplace Samples", xlab = "Value", breaks = 20)
#'
#' @export
srlaplace_truncate = function(xl, xr){
  
  res <- truncate_error_checking(xl, xr, pbgrids$srlaplace)
  xl <- res$xl; xr <- res$xr
  
  
  Upper_cumsum = .Call(C_laplace_trunc_nav, xl, xr)
  
  stopifnot(
    "xl is has a CDF close to 1" = (Upper_cumsum[1] != 1),
    "xr is has a CDF close to 0" = (Upper_cumsum[2] != 0)
  )
  
  function_string <- paste0("function(n) { .Call(C_laplace_trunc, n, ", paste0(xl), ", ", paste0(xr), ", ", paste0(Upper_cumsum[1]),
                            ", ", paste0(Upper_cumsum[2]), ", ", paste0(as.integer(Upper_cumsum[3])), ", ", paste0(as.integer(Upper_cumsum[4])), ")}")
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
}





#' @export
srlaplace_optimize = function(
    mu = 0,
    b = 1,
    steps = 4091,
    grid_range = NULL,
    theta = NULL,
    target_sample_size = 1000,
    verbose = FALSE
) {
  
  density_name <- 'srlaplace'
  
  dendata <- pbgrids[[density_name]]
  
  f_params <- c(mu, b) # F L
  
  modes <- dendata$set_modes(mu)
  
  f <- dendata$create_f(mu, b)
  
  if( identical(dendata$tails_method,"ARS") ){
    h <- function(x)
      log(f(x))
    
    h_prime <- stors_prime(modes, h)
  }else{
    cdf <- dendata$create_cdf(mu, b)
  }
  
  grid_optimizer(dendata, density_name, f, cdf, h,
                 h_prime, modes, f_params, steps,
                 grid_range, theta, target_sample_size, verbose)
  
}

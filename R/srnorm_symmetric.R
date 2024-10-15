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
#' @export
srnorm_symmetric <- function(n) {
  .Call(C_srnorm_symmetric, n)
}





#' @export
srnorm_symmetric_optimize = function(
    steps = 4091,
    grid_range = NULL,
    theta = NULL,
    target_sample_size = 1000,
    verbose = FALSE
) {
  
  mu <- 0; sd <- 1
  
  density_name <- 'srnorm_symmetric'
  
  dendata <- pbgrids[[density_name]]
  
  f_params <- c(mu, sd) # F L
  
  modes <- mu
  
  f <- dendata$f
  
  h <- dendata$h
  
  h_prime <- dendata$h_prime
  
  grid_optimizer(dendata, density_name, f, cdf, h,
                 h_prime, modes, f_params, steps,
                 grid_range, theta, target_sample_size, verbose)
  
}


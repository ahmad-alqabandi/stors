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
#' @export
srgamma <- function(n) {
  .Call(C_srgamma_custom, n)
}



#' @export
srgamma_optimize = function(shape = 1,
                            rate = 1,
                            xl = NULL,
                            xr = NULL,
                            scale = 1 / rate,
                            steps = 4091,
                            grid_range = NULL,
                            theta = NULL,
                            target_sample_size = 1000,
                            verbose = FALSE) {
  
  dist_name <- 'srgamma'
  
  symmetric <- NULL
  
  dendata <- pbgrids[[dist_name]]
  
  cnum <- dendata$Cnum + 1
  
  grid_type = "custom"
  
  f_params <- list(shape = shape, scale = scale) # F L
  
  modes <- dendata$set_modes(shape, scale)
  
  f <- dendata$create_f(shape = shape, scale = scale)
  
  grid_optimizer(
    dendata,
    dist_name,
    xl,
    xr,
    f,
    modes,
    f_params,
    steps,
    grid_range,
    theta,
    target_sample_size,
    grid_type,
    symmetric,
    cnum,
    verbose
  )
  
}

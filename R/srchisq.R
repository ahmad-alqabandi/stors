#' Sampling from Normal Distribution
#' @rdname srchisq
#' @order 1
#'
#' @description
#' Sampling from the Normal distribution using Stors.
#'
#' @details
#' The function \code{srchisq()} is used for sampling from a standard Normal distribution (mean = 0, standard deviation = 1).
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
#' \code{srchisq()} returns a sample of size \code{n} from a standard normal distribution.
#'
#' @examples
#' # Generating Samples from a Standard Normal Distribution
#' # This example illustrates how to generate 10 samples from a standard normal distribution.
#' # It first optimizes the grid for sampling using \code{grid_optimizer},
#' # and then generates samples using \code{srchisq}.
#'
#' # Optimize the grid for the standard normal distribution
#' grid_optimizer("srchisq")
#'
#' # Generate and print 10 samples from the standard normal distribution
#' samples <- srchisq(10)
#' print(samples)
#'
#'
#' @export
srchisq <- function(n = 1, x = NULL) {
  .Call(C_srchisq_custom_check, n,x)
}




#' @export
srchisq_optimize = function(
  df = NULL,
  xl = NULL,
  xr = NULL,
  steps = 4091,
  grid_range = NULL,
  theta = NULL,
  target_sample_size = 1000,
  verbose = FALSE
) {
  
  if(is.null(df)) df <- 2
  
  dist_name <- 'srchisq'
  
  symmetric <- NULL
  
  dendata <- pbgrids[[dist_name]]
  
  cnum <- dendata$Cnum + 1
  
  grid_type = "custom"
  
  if (df == 1) {
    message("Grid building is not available for df = 1. You can square the result of srnorm() to sample from this distribution.")
    return()
  }
  
  f_params <- list(df = df)
  
  modes <- dendata$set_modes(df)

  f <- dendata$create_f(df)
  
  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type,
                 symmetric,
                 cnum, verbose)
}

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
srexp <- function(n, rate = 1) {
  .Call(C_srexp_scaled, n, c(rate))
}

#' Sampling from Custom Exponential Distribution
#' @rdname srexp
#' @order 2
#' @export
srexp_custom <- function(n) {
  .Call(C_srexp_custom, n)
}



#' @export
srexp_optimize = function(
    rate = 1,
    xl = NULL,
    xr = NULL,
    steps = 4091,
    grid_range = NULL,
    theta = NULL,
    target_sample_size = 1000,
    verbose = FALSE
    ) {
  
  dist_name <- 'srexp'
  symmetric <- NULL
  
  dendata <- pbgrids[[dist_name]]
  
  if(dendata$scalable){
    if(identical(list(rate = rate), dendata$std_params))
    {
      cnum <- dendata$Cnum
      grid_type = "scaled"
    } else {
      cnum <- dendata$Cnum + 1
      grid_type = "custom"
    }
  }else{
    cnum <- dendata$Cnum + 1
    grid_type = "custom"
    
  }
  
  f_params <- list(rate = rate)
  
  modes <- 0
  
  f <- dendata$create_f(rate)
  
  
  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size, grid_type,
                 symmetric, cnum, verbose)
  
}

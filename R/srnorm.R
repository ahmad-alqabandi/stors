


#' Sampling from Normal distribution
#' @rdname snorm
#' @order 1
#' 
#' @description
#' Sampling from Normal distribution using stors.
#' 
#' @details
#' The function srnorm() is used to sample from a standard Normal distribution with a mean equal to zero and a standard deviation equal to 1. In contrast, srsnorm() is used to sample from any Normal distribution by accepting the mean and standard deviation as inputs. The reason for having these two separate functions is to improve performance. The Stors algorithm is so fast that even simple operations like multiplication and addition could impact its speed.
#' 
#' @param n sample size
#' 
#' @return sample of size n
#' 
#' 
#' @example
#' # the following example shows how to generate 10 samples from standard normal distribution using stors
#' 
#' grid_optimizer("srnorm")
#' 
#' srnorm(10)
#' 
#' # the following example shows how to generate 10 samples from normal distribution with mean equal 4 and standard deviation equal to 2
#'
#' srsnorm(10,4,2)
#' 
#' @export
srnorm <- function(n) {
  .Call(C_srnorm, n)
}



#' @rdname snorm
#' @order 2
#' 
#' @param n Integer sample size
#' @param mu Scalar mean
#' @param sd Scalar standard deviations
#'
#' @export
srsnorm <- function(n, mu =0, sd = 1) {
  .Call(C_srnorm, n) * sd + mu
}


#' Title
#'
#' @export
d_srnorm_upper = function(n ,xl, xr ,csl , csr, il, ir){
  
  .Call(C_srnorm_trunc, n,xl, xr, csl, csr, il, ir)
  
}


# I will stick with this function until I found a proper solution.
#' Title
#'
#' @param n 
#' @param l 
#' @param r 
#'
#' @return
#' @export
#'
#' @examples
truncsrnorm = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >=  pbgrids$srnorm$lb,
    "xr must be smaller than the density upper bound" = xr <=  pbgrids$srnorm$rb
  )
  
  Upper_cumsum = .Call(C_srnorm_trunc_nav, xl, xr)
  
  print(Upper_cumsum)
  return(
    function(n){
      d_srnorm_upper(n, xl, xr, Upper_cumsum[1], Upper_cumsum[2], as.integer(Upper_cumsum[3]), as.integer(Upper_cumsum[4]))
    }
  )
  
}

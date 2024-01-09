
#' laplace
#'
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
laplace <- function(n) {
  .Call(C_laplace, n)
}


#' Title
d_laplace_upper = function(n ,xl, xr ,csl , csr, il, ir){
  
  .Call(C_laplace_trunc, n,xl, xr, csl, csr, il, ir)
  
}


#' Title
#' @param n 
#' @param l 
#' @param r 
#' @return
#' @export
#'
#' @examples
trunclaplace = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >=  pbgrids$laplace$lb,
    "xr must be smaller than the density upper bound" = xr <= pbgrids$laplace$rb
  )
  
  Upper_cumsum = .Call(C_laplace_trunc_nav, xl, xr)
  
  print(Upper_cumsum)
  return(
    function(n){
      d_laplace_upper(n, xl, xr, Upper_cumsum[1], Upper_cumsum[2], as.integer(Upper_cumsum[3]), as.integer(Upper_cumsum[4]))
    }
  )
}



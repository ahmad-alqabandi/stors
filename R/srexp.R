#' @export
srexp <- function(n) {
  .Call(C_srexp, n)
}




#' Title
#'
d_srexp_upper = function(n ,xl, xr ,csl , csr, il, ir){
  
  .Call(C_srexp_trunc, n,xl, xr, csl, csr, il, ir)
  
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
truncsrexp = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >=  pbgrids$srexp$lb,
    "xr must be smaller than the density upper bound" = xr <= pbgrids$srexp$rb
  )
  
  Upper_cumsum = .Call(C_srexp_trunc_nav, xl, xr)
  
  print(Upper_cumsum)
  return(
    function(n){
      d_srexp_upper(n, xl, xr, Upper_cumsum[1], Upper_cumsum[2], as.integer(Upper_cumsum[3]), as.integer(Upper_cumsum[4]))
    }
  )
  
}

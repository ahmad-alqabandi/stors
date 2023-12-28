#' laplace
#'
trunclaplace_inter = function(n ,l, r){
  
  .Call(C_laplace_trunc, n, l, r)
  
}



# I will stick with this function until I found a proper selution.
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
trunclaplace = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >  pbgrids$srnorm$lb,
    "xr must be smaller than the density upper bound" = xr < pbgrids$srnorm$rb
  )
  
  Upper_cumsum = .Call(C_laplace_trunc_nav, xl, xr)
  
  return(
    function(n){
      trunclaplace_inter(n, Upper_cumsum[1], Upper_cumsum[2])
    }
  )
  
}






# Using paste0 which makes it slower, what is happening here is when the parent function get called it fetches the values of CDF(x_l) and CDF(x_r) from memmory then apply function paste0() which make it slower anstede.
#' @export
truncsrnorm_paste_0 = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >  pbgrids$srnorm$lb,
    "xr must be smaller than the density upper bound" = xr < pbgrids$srnorm$rb
  )
  
  Upper_cumsum = .Call(C_srnorm_trunc_nav, xl, xr)

  return(
    function(n){
      truncsrnorm_inter(n, paste0(Upper_cumsum[1]), paste0(Upper_cumsum[2]))
      # truncsrnorm_inter(n, 0.02325644, 0.99846791)
      
    }
  )

}








# Hard coded just as a standard measure
#' @export
truncsrnorm_hard_code = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >  pbgrids$srnorm$lb,
    "xr must be smaller than the density upper bound" = xr < pbgrids$srnorm$rb
  )
  
  Upper_cumsum = .Call(C_srnorm_trunc_nav, xl, xr)

  return(
    function(n){
      # truncsrnorm_inter(n, Upper_cumsum[1], Upper_cumsum[2])
      truncsrnorm_inter(n, 0.02325644, 0.99846791)
      
    }
  )
  
}




# Here we parse the function to hard code the values of CDF(x_l) and CDF(x_r).
#' @export
truncsrnorm_trick = function(xl, xr){
  
  stopifnot(
    "xl must be smaller that xr" = xl < xr,
    "xl must be greater than the density lower bound" = xl >  pbgrids$srnorm$lb,
    "xr must be smaller than the density upper bound" = xr < pbgrids$srnorm$rb
  )
  
  Upper_cumsum = .Call(C_srnorm_trunc_nav, xl, xr)

  function_string <- paste0("function(n){  truncsrnorm_inter(n, ",paste0( Upper_cumsum[1]),", ",paste0( Upper_cumsum[1]),") }")
  
  function_expression <- parse(text = function_string)
  
  my_function <- eval(function_expression)
  
  
  return(as.function(my_function))
  
}






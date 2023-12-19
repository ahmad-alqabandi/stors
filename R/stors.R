#' Sampling Function for Users' Grid
#'
#' @description
#' This function generates a sampling function based on a grid created and optimized by the user using the `build_grid()` function. The resulting sampling function can then be used to produce samples.
#'
#' @param grid The sampling grid created by the user.
#'
#' @return 
#' Returns a sampling function that can be used to generate samples from the input `grid`.
#' 
#' @import digest digest
#' @export
stors <- function(grid) {
  
  force(grid)
  
  is_valid_grid(grid)
  
  if(digest(grid) %in% stors_env$user_cached_grids$Id ){
    
        Cnum = subset(stors_env$user_cached_grids, Id == digest(grid))$Cnum
    
  } else{
    
    n = nrow(stors_env$user_cached_grids) + 1
    
    stors_env$user_cached_grids[n,]$Id = digest(grid)
    
    Cnum = stors_env$user_cached_grids[n,]$Cnum = stors_env$grids$biultin$builtin_num + n - 1
    
    cash_grid_c(Cnum, grid)
    
  }
  
  dens_func <- eval(parse(text = grid$dens_func))
  
  rfunc_env <- new.env()
  
  rm(grid)
  
  function(n) {
    .Call(
      C_stors,
      n,
      eval(expression(Cnum)),
      dens_func,
      rfunc_env
    )
  }
}


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



#' Sampling from Laplace distribution
#'
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
srlaplace <- function(n) {
  .Call(C_slaplace, n)
}


#' Sampling from Laplace distribution
#'
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
rlaplace_c <- function(n) {
  .Call(C_rLaplace_c, n)
}

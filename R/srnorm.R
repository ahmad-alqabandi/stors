
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
#' \deqn{ f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mean)^2}{2\sigma^2}} }
#' 
#' Where \eqn{\mean} is the mean and \eqn{\sigma} is the standard deviation.
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
#'
#' @export
srnorm <- function(n, mean = 0, sd = 1) {
  # C_scaled_sampling_fun <- get("C_scaled_sampling_fun", envir = srnorm_fun_env)
  .Call(C_srnorm_scaled_check, n, c(mean, sd))
}

#' @export
srnorm_custom <- function(n) {
  .Call(C_srnorm_custom_check, n)
}

#' @rdname srnorm
#' @order 2
#'
#' @param xl Lower bound for truncation.
#' @param xr Upper bound for truncation.
#' 
#' 
#' @details
#' \code{srnorm_truncate()}, this function allows sampling from a standard normal distribution that is truncated within specified bounds.
#'  It is particularly useful when the area of interest in a normal distribution is limited to a specific range.
#'  The function first validates the truncation bounds to ensure they are within the allowable range of the distribution and then creates a tailored sampling function based on these bounds.
#' 
#' @return
#' \code{srnorm_truncate()} returns a function that, when called with a sample size \code{n}, generates \code{n} samples from a normal distribution truncated between \code{xl} and \code{xr}.
#'
#' @examples
#' # Generating Samples from a Truncated Standard Normal Distribution
#' # This example demonstrates how to generate 100,
#' # samples from a standard normal distribution truncated in the range [-2, 2].
#'
#' # Create the truncated sampling function
#' norm_trunc <- srnorm_truncate(xl = -2, xr = 2)
#'
#' # Generate 100 samples
#' sample <- norm_trunc(100)
#'
#' # Plot a histogram of the samples
#' hist(sample, main = "Histogram of Truncated Normal Samples", xlab = "Value", breaks = 20)
#' 
#' @export
srnorm_truncate <- function(xl = -Inf, xr = Inf, mean = 0, sd = 1){
  
  dist_name <- 'srnorm'
  
  dendata <- pbgrids[[dist_name]]
  
  cnum_scalable <- dendata$Cnum
  cnum_custom <- dendata$Cnum + 1
  choosen_grid_num <- NULL
  
  scalable_info <-  cached_grid_info(cnum_scalable)
  custom_info <-  stors:::cached_grid_info(cnum_custom)
  
  res <- truncate_error_checking(xl, xr, dendata)
  xl <- res$xl; xr <- res$xr
  
  if( !is.null(custom_info) && all(custom_info[-1] == c(mean, sd)) ){
    choosen_grid_num <- cnum_custom
    Upper_cumsum <- .Call(C_srnorm_trunc_nav, xl, xr, choosen_grid_num)
    
  }else if( !is.null(scalable_info)){
    choosen_grid_num <- cnum_scalable
    Upper_cumsum <- .Call(C_srnorm_trunc_nav, xl, xr, choosen_grid_num)
    
  } else{
    .Call(C_grid_error,0,0)
    return(NULL)
    
  }
  

  stopifnot(
    "xl is has a CDF close to 1" = (Upper_cumsum[1] != 1),
    "xr is has a CDF close to 0" = (Upper_cumsum[2] != 0)
  )
  
  function_string <- paste0("function(n) { .Call(C_srnorm_trunc, n, ",
                            paste0(xl), ", ", paste0(xr), ", ", paste0(Upper_cumsum[1]),
                            ", ", paste0(Upper_cumsum[2]), ", ", paste0(as.integer(Upper_cumsum[3])), ", ",
                            paste0(as.integer(Upper_cumsum[4])),", ", 
                            paste0(as.integer(choosen_grid_num)), ")}")
  
  function_expression <- parse(text = function_string)
  sampling_function <- eval(function_expression)
  
  return(sampling_function)
}





#' @export
srnorm_optimize = function(
  mean = 0,
  sd = 1,
  xl = NULL,
  xr = NULL,
  steps = 4091,
  grid_range = NULL,
  theta = NULL,
  target_sample_size = 1000,
  verbose = FALSE,
  symmetric = NULL) {
  
  dist_name <- 'srnorm'
  
  dendata <- pbgrids[[dist_name]]

  if(dendata$scalable){
    if(identical(list(mean = mean, sd = sd), dendata$std_params))
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
  
  f_params <- list(mean = mean, sd = sd)
  
  modes <- dendata$set_modes(mean)
  
  f <- dendata$create_f(mean, sd)

  

  check_grid_optimization_criteria(symmetric, cnum, dendata)
  
  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type,
                 symmetric,
                 cnum, verbose)
  
}

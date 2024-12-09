
#' Sampling from the Normal Distribution
#' @rdname srnorm
#' @order 1
#'
#' @description
#' The \code{srnorm()} function generates random samples from a Normal distribution using the STORS algorithm. 
#' It employs an optimized proposal distribution around the mode and Adaptive Rejection Sampling (ARS) for the tails.
#'
#' @details
#' By default, \code{srnorm()} samples from a standard Normal distribution (\code{mean = 0}, \code{sd = 1}). 
#' The proposal distribution is pre-optimized at package load time using \code{srnorm_optimize()} with 
#' \code{steps = 4091}, creating a scalable grid centered around the mode. 
#'
#' If \code{srnorm()} is called with custom \code{mean} or \code{sd} parameters, the samples are generated 
#' from the standard Normal distribution and scaled accordingly.
#'
#' When \code{srnorm_optimize()} is explicitly called:
#' - A grid is created and cached. If no parameters are provided, a standard grid is created (\code{mean = 0}, \code{sd = 1}).
#' - Providing \code{mean} or \code{sd} creates a custom grid, which is cached for use with \code{srnorm_custom()}.
#' - The optimization process can be controlled via parameters such as \code{steps}, \code{grid_range}, or 
#'   \code{theta}. If no parameters are provided, the grid is optimized via brute force based on the 
#'   \code{target_sample_size}.
#'
#' For efficiency, \code{srnorm_optimize()} also supports one-tailed grids via the \code{symmetric} argument. 
#' This reduces memory usage for symmetrical distributions, though sampling may be slower for large samples 
#' due to the additional computation required for symmetry checks.
#'
#' The Normal distribution has the density:
#' 
#' \deqn{ f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mu)^2}{2\sigma^2}} }
#' 
#' where \eqn{\mu} is the mean and \eqn{\sigma} is the standard deviation.
#'
#' @param n Integer. Number of samples to draw.
#' @param mean Numeric. Mean of the target Normal distribution.
#' @param sd Numeric. Standard deviation of the target Normal distribution.
#' @param x Numeric vector of length \eqn{n}. If provided, samples are stored in this vector to avoid memory 
#' allocation, improving performance for repeated sampling of the same size.
#'
#' @return
#' A numeric vector of length \code{n} containing samples from the Normal distribution with the specified 
#' \code{mean} and \code{sd}.
#'
#' @examples
#' # Generate 10 samples from the standard Normal distribution
#' samples <- srnorm(10)
#' print(samples)
#'
#' # Generate 10 samples using a pre-allocated vector
#' x <- numeric(10)
#' srnorm(10, x = x)
#' print(x)
#'
#' # Generate 10 samples from a Normal distribution with mean = 2 and sd = 3
#' samples <- srnorm(10, mean = 2, sd = 3)
#' print(samples)
#'
#' @export
srnorm <- function(n = 1, mean = 0, sd = 1, x = NULL) {
  .Call(C_srnorm_scaled_check, n, c(mean, sd), x)
}


#' Sampling from Custom Normal Distribution
#' @rdname srnorm
#' @order 2
#' 
#' @export
srnorm_custom <- function(n = 1, x = NULL) {
  .Call(C_srnorm_custom_check, n, x)
}


#' Optimizing Normal Distribution Grid
#' @rdname srnorm
#' @order 3
#' 
#' @param xl Scalar lower bounds.
#' @param xr Scalar upper bounds.
#' @param steps Optional Scalar integer indicating the number of steps in the proposal distribution.
#' @param grid_range Optional Vector of two elements specifying the start and end points for constructing steps along the x-axis.
#' @param theta Optional Scalar defining the pre-acceptance threshold,
#'  dictating when the proposal steps constructing break based on the probability of pre-acceptance.
#' @param target_sample_size Scalar integer indicating the target sample size. The grid optimization process will take this number into account.
#' @param verbose Boolean if set to True, a table detailing the optimization areas and steps will be displayed during grid optimization.
#' @param symmetric Scalar can be set to distribution's mode to sample from one tiled grid.
#' 
#' @export
srnorm_optimize = function(
  mean = NULL,
  sd = NULL,
  xl = -Inf,
  xr = Inf,
  steps = 4091,
  grid_range = NULL,
  theta = NULL,
  target_sample_size = 1000,
  verbose = FALSE,
  symmetric = NULL) {
  
  dist_name <- 'srnorm'
  
  dendata <- pbgrids[[dist_name]]

  f_params <- list(mean = mean, sd = sd)
  
  if(dendata$scalable){
    
    isnull <- sapply(f_params, is.null)
    
    if(all(isnull)){
      cnum <- dendata$Cnum
      grid_type = "scaled"
    }else{
      cnum <- dendata$Cnum + 1
      grid_type = "custom"
    }
    
    f_params <- ifelse(isnull, dendata$std_params, f_params)

  }else{
    cnum <- dendata$Cnum + 1
    grid_type = "custom"
  }
  
  modes <- dendata$set_modes(f_params$mean)
  
  f <- dendata$create_f(f_params$mean, f_params$sd)
  

  check_grid_optimization_criteria(symmetric, cnum, dendata)
  
  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type, symmetric, cnum, verbose)
  
}
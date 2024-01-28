#' Sampling Function for Users' Grid
#'
#' @description
#' This function generates a sampling function based on a grid created by the user using the `build_grid()` function. 
#' The resulting sampling function can then be used to produce samples, including those from truncated distributions.
#'
#' @param grid The sampling grid created using the `build_grid()` function.
#' @param xl The lower bound for truncation. Default is the left bound of the grid.
#' @param xr The upper bound for truncation. Default is the right bound of the grid.
#'
#' @return 
#' Returns a function that can be used to generate samples from the specified `grid`. If `xl` and `xr` are provided, 
#' the samples are drawn from the truncated distribution within these bounds.
#'
#' @details
#' After a user creates a proposal grid for their desired sampling function using `build_grid()`,
#' this grid must be passed to `stors()` to create a sampling function for the target distribution. 
#' `stors()` first checks whether the grid was indeed created using `build_grid()`. If the user has altered 
#' or modified the grid returned from `build_grid()`, `stors()` will reject the altered grid; therefore, 
#' no changes should be made to the grid after its creation. Once the grid is accepted by `stors()`, it is 
#' cached in memory, allowing fast access to grid data for the compiled C code and reducing memory access latency. 
#' Subsequently, `stors()` returns a function that can be utilized to generate samples from the target distribution,
#' optionally truncated between `xl` and `xr`. The truncation in `stors` is achieved by truncating the proposal grid itself,
#' ensuring that samples are always produced within the truncation range [xl, xr]. Therefore, if the truncation values fall within
#' the steps range of the proposal grid, no integration of the density function is required.
#'
#' @examples
#' 
#' # Example 1
#' # To sample from a standard normal distribution \( f(x) \sim \mathcal{N}(0,1) \), first build the proposal grid using `build_grid()`
#'
#' modes_norm = 0
#' f_norm <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x^2) }
#' h_norm <- function(x) { log(f_norm(x)) }
#' h_prime_norm <- function(x) { -x }
#' normal_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_norm, f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)
#'
#' # Generate samples from the standard normal distribution
#' sample_normal <- stors(normal_grid)
#' hist(sample_normal(100), main = "Normal Distribution Samples")
#'
#' # To sample from a truncated standard normal distribution between -1 and 1
#' sample_truncated_normal <- stors(normal_grid, xl = -1, xr = 1)
#' hist(sample_truncated_normal(100), main = "Truncated Normal Distribution Samples")
#'
#'
#' # Example 2
#' # Let's consider a bimodal distribution composed of two normal distributions:
#' #The first normal distribution N(0,1) with weight p = 0.3,
#' #and the second normal distribution N(4,1) with weight q = 0.7.
#' 
#' f_bimodal <- function(x) {
#'  0.3 * (1 / sqrt(2 * pi) * exp(-0.5 * (x - 0)^2)) +
#'  0.7 * (1 / sqrt(2 * pi) * exp(-0.5 * (x - 4)^2))
#'}
#'
#' # Define the modes of the bimodal distribution
#'    modes_bimodal <- c(0.00316841, 3.99942)
#'
#' # Build the proposal grid for the bimodal distribution
#' bimodal_grid = build_grid(f = f_bimodal, modes = modes_bimodal, lb = -Inf, rb = Inf, steps = 1000)
#' 
#' # Create the sampling function using `stors()`
#' sample_bimodal <- stors(bimodal_grid)
#' 
#' # Generate and plot samples from the bimodal distribution
#' bimodal_samples <- sample_bimodal(1000)
#' hist(bimodal_samples, breaks = 30, main = "Bimodal Distribution Samples")
#' 
#' # Create the truncated sampling function using `stors()` with truncation bounds [0.5, 6]
#' sample_truncated_bimodal <- stors(bimodal_grid, xl = 0.5, xr = 6)
#'
#' # Generate and plot samples from the truncated bimodal distribution
#' truncated_bimodal_samples <- sample_truncated_bimodal(1000)
#' hist(truncated_bimodal_samples, breaks = 30, main = "Truncated Bimodal Distribution Samples")
#' 
#' @import digest digest
#' @export
stors <- function(grid, xl = grid$grid_bounds[1], xr = grid$grid_bounds[2]) {
  

  force(grid)
  
  is_valid_grid(grid)

  Cnum <- cache_stors_grid(grid)
  
  dens_func <- eval(parse(text = grid$dens_func))
  
  rfunc_env <- new.env()
  
  
  if( xl != grid$grid_bounds[1] || xr != grid$grid_bounds[2]){

      stopifnot(
        "xl must be smaller that xr" = xl < xr,
        "xl must be greater than or equal the density lower bound" = xl >  grid$grid_bounds[1],
        "xr must be smaller than or equal the density upper bound" = xr <  grid$grid_bounds[2]
      )

    Upper_cumsum = .Call(C_stors_trunc_nav,Cnum, xl, xr)

    function_string <- paste0("function(n) { .Call(C_stors_trunc, n, ",paste0(Cnum),", ",paste0(xl),", ",paste0(xr),", ",paste0(Upper_cumsum[1]),",", paste0(Upper_cumsum[2]),", ", paste0(as.integer(Upper_cumsum[3])),", ", paste0(as.integer(Upper_cumsum[4]))," , dens_func, rfunc_env) }")
    
  }else{

    function_string <- paste0("function(n) { .Call(C_stors, n, ",paste0(Cnum),", dens_func, rfunc_env) }" )
  }

  function_expression <- parse(text = function_string)
  
  sampling_function <- eval(function_expression)
  
  rm(grid)
  
  return(sampling_function)
}


#' #' @import digest digest
#' trunc_stors <- function(grid, xl, xr) {
#' 
#' 
#'   force(grid)
#' 
#'   is_valid_grid(grid)
#' 
#'   stopifnot(
#'     "xl must be smaller that xr" = xl < xr,
#'     "xl must be greater than or equal the density lower bound" = xl >  grid$grid_bounds[1],
#'     "xr must be smaller than or equal the density upper bound" = xr <  grid$grid_bounds[2]
#'   )
#' 
#'   Cnum <- cache_stors_grid(grid)
#' 
#'   dens_func <- eval(parse(text = grid$dens_func))
#' 
#'   rfunc_env <- new.env()
#' 
#'   rm(grid)
#' 
#'   Upper_cumsum = .Call(C_stors_trunc_nav,Cnum, xl, xr)
#' 
#'   print(Upper_cumsum)
#' 
#'   function_string <- paste0("function(n) { .Call(C_stors_trunc, n, ",paste0(Cnum),", ",paste0(xl),", ",paste0(xr),", ",paste0(Upper_cumsum[1]),",", paste0(Upper_cumsum[2]),", ", paste0(as.integer(Upper_cumsum[3])),", ", paste0(as.integer(Upper_cumsum[4]))," , dens_func, rfunc_env) }")
#'   function_expression <- parse(text = function_string)
#'   sampling_function <- eval(function_expression)
#' 
#'   return(sampling_function)
#' 
#' 
#' }



cache_stors_grid <- function(grid){
  
  if(digest(grid) %in% stors_env$user_cached_grids$Id ){
    
    Cnum = subset(stors_env$user_cached_grids, Id == digest(grid))$Cnum
    
  } else{
    
    n = nrow(stors_env$user_cached_grids) + 5
    # 5 is padding
    
    stors_env$user_cached_grids[n,]$Id = digest(grid)
    
    Cnum <- stors_env$user_cached_grids[n,]$Cnum <- stors_env$grids$biultin$builtin_num + n
    
    cache_grid_c(Cnum, grid)
    
  }
  
  return(Cnum)
  
}

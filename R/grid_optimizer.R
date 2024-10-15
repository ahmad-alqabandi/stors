
#' Optimize Built-in Grid
#'
#' @description
#' This function optimizes and stores the sampling grid for STROS built-in distributions.
#'
#' @param steps Optional Scalar integer indicating the number of steps in the proposal distribution.
#' @param grid_range Optional Vector of two elements specifying the start and end points for constructing steps along the x-axis.
#' @param theta Optional Scalar defining the pre-acceptance threshold, dictating when the proposal steps constructing break based on the probability of pre-acceptance.
#' @param target_sample_size Scalar integer indicating the target sample size. The grid optimization process will take this number into account.
#' @param verbose Boolean if set to True, a table detailing the optimization areas and steps will be displayed during grid optimization.
#'
#' @details
#' Before sampling from built-in distributions in the STROS package,
#'  it is necessary to first optimize the sampling grid for the specific distribution using \code{grid_optimizer()}.
#'  This function calculates and stores an optimized grid based on the given parameters or defaults.
#'  The optimization primarily focuses on balancing memory usage and computational efficiency,
#'  ensuring that the grid is both space-efficient and capable of facilitating fast sampling.
#'
#' For optimal grid configuration, it is advised that users refrain from manually setting the \code{Theta}, \code{steps}, or \code{grid_range} parameters, as \code{grid_optimizer()} automatically determines the best settings.
#'  However, users can provide their own values for these parameters to customize the grid according to specific requirements. For more detailed information about these arguments and optimization process, please refer to the manual page of \code{build_final_grid()}.
#'
#'
#' @examples
#' # Optimize grid for sampling from the normal distribution
#' library(stors)
#' grid_optimizer("srnorm")
#' # Generate 10 samples from the optimized standard normal distribution
#' srnorm(10)
#'
#' @return
#' This function generates and stores and return an optimized grid for the specified built-in distribution in R's internal data directory.
#'\item{"grid_data"}{A data frame including the created steps information, such as \code{x} (the beginning of each step on the x-axis), \code{s_upper} (the step height on the y-axis), \code{p_a} (pre-acceptance probability for each step), and \code{s_upper_lower} (a vector used to re-scale the uniform random number when the sample is accepted).}
#'\item{"areas"}{A vector containing the areas under the left tail bound, the steps in the middle, and the right tail bound.}
#'\item{"steps_number"}{A scalar representing the number of steps in the proposal.}
#'\item{"sampling_probabilities"}{A vector containing the areas under the left tail and the combined area of the left tail and middle steps.}
#'\item{"unif_scaler"}{A scalar representing the inverse probability of sampling from the step part of the proposal, \eqn{\frac{1}{p(lb < x < rb)}}. Similar to 's_upper_lower' in the 'grid_data' data frame, this value is used to scale the uniform random value when sampling from the steps part of the proposal.}
#'\item{"lt_properties"}{A vector including 5 values used when sampling under the proposal's left tail using the ARS (Adaptive Rejection Sampling) method.}
#'\item{"rt_properties"}{A vector including 6 values used when sampling under the proposal's right tail using the ARS method.}
#'\item{"alpha"}{A scalar representing the uniform step area.}
#'\item{"tails_method"}{A string representing the tails sampling method, either 'ARS' for Adaptive Rejection Sampling or 'IT' for Inverse Transform.}
#'\item{"grid_bounds"}{A vector including the left and right bounds of the target density.}
#'\item{"dens_func"}{The function passed by the user for the target density \code{f}.}
#' 
#'
#' @export
grid_optimizer <- function(dendata,density_name,
                           f,
                           cdf = NULL,
                           h = NULL,
                           h_prime = NULL,
                           modes,
                           f_params = NULL,
                           steps = NULL,
                           grid_range = NULL,
                           theta = NULL,
                           target_sample_size = NULL,
                           verbose = FALSE) {
  # density_name <- match.arg(density_name)
  # 
  # dendata <- pbgrids[[density_name]]
  free_cache_cnum_c(dendata$Cnum)
  
  grid_param <- list(
    target = list(
      density = f,
      log_density = NULL,
      log_density_prime = NULL,
      Cumulitive_density = NULL,
      modes = modes,
      modes_count = length(modes),
      between_minima = NULL,
      right_bound = dendata$rb,
      left_bound = dendata$lb
    ),
    proposal = list(
      grid_range = grid_range,
      tails_method = dendata$tails_method,
      steps = steps,
      optimal_step_area = NULL,
      left_steps_proportion = NULL,
      right_steps_proportion = NULL,
      pre_acceptance_threshold = theta,
      target_sample_size = target_sample_size
    ),
    built_in = TRUE,
    cnum = dendata$Cnum,
    symmetric = dendata$symmetric,
    c_function_name = density_name,
    verbose = verbose
  )
  
  if (grid_param$proposal$tails_method == "IT") {
    grid_param$target$Cumulitive_density = cdf
  } else{
    grid_param$target$log_density = h
    grid_param$target$log_density_prime = h_prime
  }
  
  grid_param <- grid_error_checking_and_preparation(grid_param)
  
  optimal_grid_params = find_optimal_grid(grid_param)
  opt_grid <- build_final_grid(gp = optimal_grid_params)
  opt_grid$density_parameters <- f_params
  
  cache_grid_c(dendata$Cnum, opt_grid)
  save_builtin_grid(dendata$Cnum, opt_grid)
  stors_env$grids$builtin[[density_name]]$opt <- TRUE
  
  #opt_grid$dens_func <- deparse(f)
  opt_grid$dens_func <- f
  
  
  class(opt_grid) <- "grid"
  
  return(opt_grid)
  
}


#' @importFrom microbenchmark microbenchmark
find_optimal_grid <- function(gp) {
  
  target_sample_size <- gp$proposal$target_sample_size
  theta <- gp$proposal$pre_acceptance_threshold
  grid_range <- gp$proposal$grid_range
  lb <- gp$target$left_bound
  rb <- gp$target$right_bound
  modes <- gp$target$modes
  f <- gp$target$density
  lstpsp <- gp$proposal$left_steps_proportion
  rstpsp <- gp$proposal$right_steps_proportion
  mode_n <- gp$target$modes_count
  steps <- gp$proposal$steps
  c_function_name <- gp$c_function_name
  verbose <- gp$verbose
  
  times = ceiling(opt_times / target_sample_size)
  
  if ((theta == 0 &&
       identical(grid_range, c(lb, rb))) || !is.null(steps))
  {
    max_left_stps = find_left_steps(
      gp = gp,
      area = 0.001,
      mode_i = 1,
      steps_lim = Inf
    )$steps
    max_right_stps = find_right_steps(
      gp = gp,
      area = 0.001,
      mode_i = mode_n,
      steps_lim = Inf
    )$steps
    lstpsp = max_left_stps / (max_left_stps + max_right_stps)
    rstpsp = 1 - lstpsp
    if (!is.null(steps)) {
      gp$proposal$optimal_step_area = 1 / steps
      gp$proposal$left_steps_proportion = lstpsp
      gp$proposal$right_steps_proportion = rstpsp
      return(gp)
    }
  }
  
    performance = data.frame(area = numeric(),
                             time = numeric(),
                             steps = numeric())
    
    if (gp$built_in) {
      cnum <- gp$cnum
      density_fun <- get(gp$c_function_name, mode = "function")
    } else{
      cnum <- 0
      rfunc_env <- new.env()
      density_fun <- function(n) {
        .Call(C_stors, n, cnum, f, rfunc_env)
      }
    }
 
    for (i in (1:length(opt_areas))) {
      area = opt_areas[i]
      
      if (is.null(rstpsp)) {
        step = Inf
      } else{
        step = opt_steps[i]
      }
      
      if (verbose) {
        cat(
          "\n=====================================\n",
          sprintf("Step: %10s | Area: %10.9f", step, area),
          "\n-------------------------------------\n",
          "       Area       |    Best Sim Time\n",
          "-------------------------------------\n"
        )
      }
      
      area_seq = seq(from = area * 0.95 ,
                     to = area * 1.05,
                     length.out = opt_alpha_length)
      steps_time = double()
      
      for (j in (1:length(area_seq))) {
        grid <- build_final_grid(gp = gp, opt_area = area_seq[j])
        cache_grid_c(cnum, grid)
        gc = gc()
        suppressWarnings({
          cost <- microbenchmark::microbenchmark(st = density_fun(target_sample_size), times = times)
        })
        free_cache_cnum_c(cnum)
        steps_time = append(steps_time, stats::median(cost$time[cost$expr == "st"]))
        if (verbose) {
          cat(sprintf("--- %10.9f ---   --- %10.2f ---\n", area_seq[j], steps_time[j]))
        }
      }
      
      if (i != 1 &&
          min(steps_time) >= min(performance$time, na.rm = TRUE)) {
        if (verbose) {
          cat("\n=====================================\n")
          cat("\nPerformance Data:\n")
          cat("     Area       |     Time     |   Steps\n")
          cat("-----------------------------------------\n")
          for (k in 1:nrow(performance)) {
            cat(
              sprintf(
                "%13.10f | %10.2f | %7d\n",
                performance$area[k],
                performance$time[k],
                performance$steps[k]
              )
            )
          }
          cat("=====================================\n\n")
        }
        break
      }
      
      
      min_ind = which(steps_time == min(steps_time, na.rm = TRUE))[1]
      performance[nrow(performance) + 1, ] = c(area_seq[min_ind], steps_time[min_ind], step)
      
    }
    
    opt_performance = performance[which(performance$time == min(performance$time, na.rm = TRUE))[1], ]
    gp$proposal$optimal_step_area = opt_performance$area
    gp$proposal$steps = opt_performance$steps
    gp$proposal$left_steps_proportion = lstpsp
    gp$proposal$right_steps_proportion = rstpsp
    
    return(gp)
    
}


opt_alpha_length = 3
opt_times = 10000
opt_cache_sizes = c(4 , 8 , 16 , 32, 64, 128, 256, 512, 1024)
opt_df_var = 32 # 4 doubles
opt_list_var = 160 # 20 doubles
opt_steps = round(((opt_cache_sizes * 1024) - opt_list_var) / opt_df_var)
opt_areas = 1 / opt_steps

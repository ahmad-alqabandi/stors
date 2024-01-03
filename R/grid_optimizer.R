


#' Optimize Built-in Grid
#'
#' @description
#' This function optimizes and stores the sampling grid for STROS built-in distributions.
#'
#' @param density_name string specifying the name of the distribution.
#'
#' @details
#' Before users can sample from a built-in distribution in the STROS package, they must first optimize the sampling grid using this function.
#'
#' @seealso [grid_optimizer()] is used to optimize and cache grid to sample from multible distrubution such as [srnorm()] and [srgamma()]
#'
#' @examples
#' # optimize grid to sample 10 values from normal distribution N(0,1) 
#' library(stors)
#' grid_optimizer("srnorm")
#' srnorm(10)
#'
#' @return
#' This function generates and stores an optimized grid for the specified built-in distribution in R's internal data directory.
#' 
#' @export
grid_optimizer <- function(density_name = stors_env$grids$biultin$names) {
  
  density_name <- match.arg(density_name)
  
  # stopifnot(" Grid already optimized for this distrubution" = !stors_env$stors_env$grids$biultin[[density_name]]$optgrids$biultin[[density_name]]$opt)
  
  dendata <- pbgrids[[density_name]]
  
  if(stors_env$grids$biultin[[density_name]]$opt){
    stors_env$grids$biultin[[density_name]]$opt = FALSE
    free_cache_cnum_c(dendata$Cnum)
    
  }
  
  opt_prob = find_optimal_grid(dendata, density_name)

  if(dendata$tails_method == "IT"){
    opt_grid <- grid_builder(a = opt_prob$area , th = 0, mode = dendata$modes, f = dendata$f, cdf = dendata$cdf, stps =opt_prob$steps , lstpsp =opt_prob$lstpsp , rstpsp= opt_prob$rstpsp)
  } else if(dendata$tails_method == "ARS"){
    opt_grid <- grid_builder(a = opt_prob$area , th = 0, mode = dendata$modes, f = dendata$f,  h = dendata$h, h_prime = dendata$h_prime, stps =opt_prob$steps , lstpsp =opt_prob$lstpsp , rstpsp= opt_prob$rstpsp)
  }
  
  cache_grid_c(dendata$Cnum, opt_grid)
  
  save_builtin_grid(dendata$Cnum, opt_grid)
  
  stors_env$grids$biultin[[density_name]]$opt <- TRUE
  
  func_to_text <- deparse(dendata$f)
  
  opt_grid$dens_func <- func_to_text
  
  class(opt_grid) <- "grid"
  
  return(opt_grid)
  
}



# {
# 
# density_name <- "srexp"
# 
# dendata <- pbgrids[[density_name]]
# 
# 
# opt_prob = find_optimal_grid(dendata, density_name)
# 
# cache = 3
# 
# opt_grid <- grid_builder(a = opt_prob$prob$area[cache] , th = 0, mode = dendata$modes, f = dendata$f,  h = dendata$h, h_prime = dendata$h_prime, stps =opt_prob$prob$steps[cache] , lstpsp =opt_prob$lstpsp , rstpsp= opt_prob$rstpsp)
# 
# cache_grid_c(dendata$Cnum, opt_grid)
# 
# 
# res = microbenchmark::microbenchmark(
#   srnorm(1000),
#   srnorm(10000),
#   srnorm(100000),
#   # srnorm(1000000),
#   times = 100
# )
# median(res$time)
# 
# free_cache_cnum_c(dendata$Cnum)
# 
# }


#' @importFrom microbenchmark microbenchmark
find_optimal_grid <- function(dendata, density_name){
  
  
  mode_n <- length(dendata$modes)
  
  left_stps = find_left_steps(dendata$lb, dendata$rb, 0.001, th=0.1, dendata$modes[1], 1, mode_n, dendata$f)$m
  right_stps = find_right_steps(dendata$lb, dendata$rb, 0.001, th=0.1, dendata$modes[mode_n], mode_n, mode_n, dendata$f)$m
  
  lstpsp = left_stps/ (left_stps+right_stps)
  rstpsp = 1 - lstpsp
  
  performance = data.frame(area = numeric(), time = numeric(), steps = numeric())
  
  for (i in (1:length(areas))) {
    
    print(cache_sizes[i])
    
    area = areas[i]
    step = steps[i]
    
    area_seq = seq(from = area * 0.9 , to= area * 1, length.out = 10)
    
    steps_time = rep(NA, length(area_seq))
    
    for (j in (1: length(area_seq))) {
      
      
      if(dendata$tails_method == "IT"){
        grid <- grid_builder(a = area_seq[j] , th = 0, mode = dendata$modes, f = dendata$f, cdf = dendata$cdf, stps =step , lstpsp =lstpsp , rstpsp= rstpsp)
      } else if(dendata$tails_method == "ARS"){
        grid <- grid_builder(a = area_seq[j] , th = 0, mode = dendata$modes, f = dendata$f,  h = dendata$h, h_prime = dendata$h_prime, stps =step , lstpsp =lstpsp , rstpsp= rstpsp)
      }
      
      density_fun <- get(density_name, mode = "function")
      
      gc = gc()
      grid$grid_data$x[1]=0.001
      cache_grid_c(dendata$Cnum, grid)
      
      cost <- microbenchmark::microbenchmark(
        st = density_fun(1000), #1000000
        times = 100
      )
      
      free_cache_cnum_c(dendata$Cnum)
      
      steps_time[j] = median(cost$time[cost$expr == "st"])
      
      
    }
    
    
    print(min(steps_time))
    
    performance[nrow(performance) + 1,] = c(area_seq[steps_time == min(steps_time)], min(steps_time),  step)

    if(min(steps_time) >  min(performance$time, na.rm = TRUE) ) break
    
    
  }
  
  return(cbind(performance[performance$time == min(performance$time, na.rm = TRUE),],lstpsp =lstpsp , rstpsp= rstpsp))
  # return(list(prob = performance,lstpsp =lstpsp , rstpsp= rstpsp))
  
}




cache_sizes = c(4 ,8 ,16 ,32, 64, 128, 256, 512, 1024)

df_var = 4

list_var = 20

steps = round( ((cache_sizes * 1024) - list_var) / df_var )

areas = 1 / steps



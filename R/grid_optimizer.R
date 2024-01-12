


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
grid_optimizer <- function(density_name = stors_env$grids$biultin$names, verbose = FALSE, target_sample_size = 1000, steps = NA) {
  
  density_name <- match.arg(density_name)
  
  stopifnot(" Grid already optimized for this distrubution" = !stors_env$grids$biultin[[density_name]]$opt)
  
  dendata <- pbgrids[[density_name]]
  
  if(is.na(steps)){
    
    opt_prob = find_optimal_grid(dendata , density_name, verbose = verbose, target_sample_size = target_sample_size)
    
  }else{
    mode_n = length(modes)
    left_stps = find_left_steps(lb = lb, rb = rb, a = 0.001, th=0.1, mode = modes[1], mode_i = 1, mode_n = mode_n, f = f)$m
    right_stps = find_right_steps(lb = lb, rb = rb, a = 0.001, th=0.1, mode = modes[mode_n], mode_i = mode_n, mode_n = mode_n, f = f)$m
    
    opt_prob = list()
    opt_prob$area = 1/steps
    opt_prob$steps = steps
    opt_prob$lstpsp =left_stps / sum(left_stps, right_stps)
    opt_prob$rstpsp= 1 - opt_prob$lstpsp
    
  }
  
  
  
  
  if(dendata$tails_method == "IT"){
    opt_grid <- grid_builder(lb = dendata$lb, rb = dendata$rb ,a = opt_prob$area , th = 0, mode = dendata$modes, f = dendata$f, cdf = dendata$cdf, stps =opt_prob$steps , lstpsp =opt_prob$lstpsp , rstpsp= opt_prob$rstpsp)
  } else if(dendata$tails_method == "ARS"){
    opt_grid <- grid_builder(lb = dendata$lb, rb = dendata$rb ,a = opt_prob$area , th = 0, mode = dendata$modes, f = dendata$f,  h = dendata$h, h_prime = dendata$h_prime, stps =opt_prob$steps , lstpsp =opt_prob$lstpsp , rstpsp= opt_prob$rstpsp)
  }
  
  # opt_grid <- grid_builder(lb = dendata$lb, rb = dendata$rb ,a = opt_prob$area , th = 0, mode = dendata$modes, f = dendata$f, cdf = dendata$cdf, stps =opt_prob$steps , lstpsp =opt_prob$lstpsp , rstpsp= opt_prob$rstpsp)
  
  cache_grid_c(dendata$Cnum, opt_grid)
  
  save_builtin_grid(dendata$Cnum, opt_grid)
  
  stors_env$grids$biultin[[density_name]]$opt <- TRUE
  
  func_to_text <- deparse(dendata$f)
  
  opt_grid$dens_func <- func_to_text
  
  class(opt_grid) <- "grid"
  
  return(opt_grid)
  
}


# 
# {
# 
# density_name <- "srnorm"
# 
# dendata <- pbgrids[[density_name]]
# 
# verbose = TRUE
# 
# modes = NA
# f = NA
# h = NA
# h_prime = NA
# cdf = NA
# 
# target_sample_size = 100
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
find_optimal_grid <- function(dendata = NA, density_name = NA , lb = NA, rb = NA, modes = NA, f = NA, h = NA, h_prime = NA, cdf = NA, verbose = FALSE, target_sample_size){
  
  times = ceiling( 100000 / target_sample_size )
  
  if(is.list(dendata)){
    lb = dendata$lb
    rb = dendata$rb
    modes = dendata$modes
    f= dendata$f
    tails_method = dendata$tails_method
    cnum = dendata$Cnum
    if(tails_method == "IT"){
      cdf = dendata$cdf
    } else{
      h = dendata$h
      h_prime = dendata$h_prime
    }
  }else{
    cnum = 0
  }
  
  # cat("\n")
  # cat("cnum = ", cnum)
  # cat("\n")
  
  mode_n <- length(modes)
  
  left_stps = find_left_steps(lb = lb, rb = rb, a = 0.001, th=0.1, mode = modes[1], mode_i = 1, mode_n = mode_n, f = f)$m
  right_stps = find_right_steps(lb = lb, rb = rb, a = 0.001, th=0.1, mode = modes[mode_n], mode_i = mode_n, mode_n = mode_n, f = f)$m
  
  lstpsp = left_stps/ ( left_stps + right_stps )
  rstpsp = 1 - lstpsp
  
  performance = data.frame(area = numeric(), time = numeric(), steps = numeric())
  
  for (i in (1:length(areas))) {
    

    area = top_areas[i]
    step = opt_steps[i]
    
    if(verbose){
      cat("\n===============================\n","\n--- steps = ", step, " --- \n\n",
          "     --- area ---   --- best sim time ---\n"
         )
    }
    
    area_seq = seq(from = area * 0.95 , to= area * 1.05, length.out = 10)
    
    #steps_time = rep(NA, length(area_seq))
    steps_time = double()
    
    for (j in (1: length(area_seq))) {
      
      
      # if(dendata$tails_method == "IT"){
      #   grid <- grid_builder( lb = dendata$lb, rb = dendata$rb ,a = area_seq[j] , th = 0, mode = dendata$modes, f = dendata$f, cdf = dendata$cdf, stps =step , lstpsp =lstpsp , rstpsp= rstpsp)
      # } else if(dendata$tails_method == "ARS"){
      #   grid <- grid_builder(lb = dendata$lb, rb = dendata$rb, a = area_seq[j] , th = 0, mode = dendata$modes, f = dendata$f,  h = dendata$h, h_prime = dendata$h_prime, stps =step , lstpsp =lstpsp , rstpsp= rstpsp)
      # }
      
      grid <- grid_builder( lb = lb, rb = rb ,a = area_seq[j] , th = 0, modes, f = f, h =  h, h_prime = h_prime , cdf = cdf, stps =step , lstpsp =lstpsp , rstpsp= rstpsp)
     
      cache_grid_c(cnum, grid)
      
      if( is.list(dendata) ){
        density_fun <- get(density_name, mode = "function")
      }else{
        rfunc_env <- new.env()
        density_fun <- function(n) { .Call(C_stors,n,0,f,rfunc_env) }
      }
      
      
      gc = gc()
      
suppressWarnings({
  
      cost <- microbenchmark::microbenchmark(
        st = density_fun(target_sample_size), 
        times = times
      )
      
})
      free_cache_cnum_c(cnum)

      steps_time = append(steps_time, median(cost$time[cost$expr == "st"]))
      
      
      if(verbose){
        cat("--- ",area_seq[j]," ---   --- ",steps_time[j]," ---\n")
      }
    }

    min_ind = which(steps_time == min(steps_time, na.rm = TRUE))[1]
    
    performance[nrow(performance) + 1,] = c(area_seq[min_ind], steps_time[min_ind],  step)

    
    if(min(steps_time) >  min(performance$time, na.rm = TRUE) && verbose){
      
      cat("\n===============================\n\n")
      print(performance)
      cat("\n===============================\n\n")
      cat("optimal grid has ", opt_steps[i-1] ," steps, whith cache size = ", cache_sizes[i-1], " Kb")
      
      break
    } 
    
    
  }
  
  return(cbind(performance[performance$time == min(performance$time, na.rm = TRUE),],lstpsp =lstpsp , rstpsp= rstpsp))
  # return(list(prob = performance,lstpsp =lstpsp , rstpsp= rstpsp))
  
}


# bm_sample_size = 100000

# bm_times = 10

opt_cache_sizes = c(4 ,8 ,16 ,32, 64, 128, 256, 512, 1024)

opt_df_var = 4

opt_list_var = 20

opt_steps = round( ((opt_cache_sizes * 1024) - opt_list_var) / opt_df_var )

top_areas = 1 / opt_steps



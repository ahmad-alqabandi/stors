build_grid_2 <- function(lb = -Inf, rb = Inf, modes, f, h = NULL,
                       h_prime = NULL, cdf = NULL, steps = NA, grid_range = c(lb, rb), theta = 0,target_sample_size = 1000, verbose = FALSE) {
  
  
  grid_prob_error_checking(steps, grid_range, theta, lb, rb, modes)
  
  if( is.null(cdf)){
    
    if(is.null(h)){
      h <- function(x){log(f(x))}
    }
    
    if(is.null(h_prime)){
      h_prime <- stors_prime(modes[1], h)
    }
    
  }
  

  
  
  if(is.na(steps)){
    
    opt_prob = find_optimal_grid(lb = lb, rb = rb, modes = modes, f = f, h =  h,
                                 h_prime = h_prime, cdf = cdf, verbose = verbose,
                                 target_sample_size = target_sample_size,
                                 theta = theta, grid_range = grid_range)
    
    opt_grid <- grid_builder(lb = lb, rb = rb , a = opt_prob$area,
                             modes, f = f, h =  h, h_prime = h_prime , cdf = cdf,
                             stps =opt_prob$steps , lstpsp =opt_prob$lstpsp , rstpsp= opt_prob$rstpsp,
                             theta = theta, grid_range = grid_range)
    
  }
  else{
    
    opt_prob = list()
    
    mode_n = length(modes)
    
    left_stps = find_left_steps(lb = lb, rb = rb, a = 0.001,
                                mode = modes[1], mode_i = 1, mode_n = mode_n, f = f,
                                theta =0.1, grid_range = grid_range)$m
    
    right_stps = find_right_steps(lb = lb, rb = rb, a = 0.001,
                                  mode = modes[mode_n], mode_i = mode_n, mode_n = mode_n, f = f,
                                  theta =0.1, grid_range = grid_range)$m
    
    lstpsp =left_stps / sum(left_stps, right_stps)
    rstpsp= 1 - lstpsp
    
    opt_grid <- grid_builder(lb = lb, rb = rb ,a = 1/steps ,
                             modes = modes, f = f, h =  h, h_prime = h_prime, cdf = cdf,
                             stps =steps , lstpsp =lstpsp , rstpsp= rstpsp,
                             theta = theta, grid_range = grid_range)
    
  }
  
  
  ftx <- deparse(f)
  
  opt_grid$dens_func <- ftx
  
  if(!is.null(cdf)){
    ftx <- deparse(cdf)
    
    opt_grid$dens_cdf <- ftx
    
    opt_grid$tails_method = "IT"
    
  }else{
    
    opt_grid$tails_method = "ARS"
    
  }
  
  
  class(opt_grid) <- "grid"
  
  
  if(!(digest(opt_grid) %in% stors_env$created_girds_Id))
    stors_env$created_girds_Id  = append(stors_env$created_girds_Id , digest(opt_grid))
  
  
  return(opt_grid)
  
}

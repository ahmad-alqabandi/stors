
#c("srnorm","slaplace","laplace")


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
  
  stopifnot(" Grid already optimized for this distrubution" = !stors_env$grids$biultin[[density_name]]$opt)
  
  dendata <- pbgrids[[density_name]]
  

  if(dendata$tails_method == "IT"){
    opt_grid <- bi_grid_builder(a = 0.0003, th = 0.6, mode = dendata$modes, f = dendata$f, cdf = dendata$cdf)
  } else if(dendata$tails_method == "ARS"){
    opt_grid <- grid_builder(a = 0.0003, th = 0.6, mode = dendata$modes, f = dendata$f,  h = dendata$h, h_prime = dendata$h_prime)
    }

  cash_grid_c(dendata$Cnum, opt_grid)

  save_builtin_grid(dendata$Cnum, opt_grid)

  stors_env$grids$biultin[[density_name]]$opt <- TRUE

  func_to_text <- deparse(dendata$f)

  opt_grid$dens_func <- func_to_text

  class(opt_grid) <- "grid"

  return(opt_grid)
  
}





bi_grid_builder <- function(lb = -Inf, rb = Inf, a, th, mode, f, cdf) {
  
  mode_n <- length(mode)
  
  final_grid <- data.frame(
    x = c(), s_upper = c(), s_lower = c(), p_a = c(),
    s_upper_lower = c()
  )
  
  grids <- list()
  area <- c(0, 0, 0)
  g_len <- c()
  
  for (mode_i in (1:mode_n)) {
    grids[[mode_i]] <- bi_find_steps(lb = -Inf, rb = Inf, a, th, mode[mode_i], mode_i, mode_n, f)
  }
  
  if (mode_n > 1) {
    for (i in (1:(mode_n - 1))) {
      if (grids[[i]]$d$x[grids[[i]]$m + 1] > grids[[i + 1]]$d$x[1]) {
        if (grids[[i]]$d$s_upper[grids[[i]]$m] > grids[[i + 1]]$d$s_upper[1]) {
          grids[[i]]$d <- head(grids[[i]]$d, -1)
          
          grids[[i]]$d$s_upper[grids[[i]]$m] <- a / (grids[[i + 1]]$d$x[1] - grids[[i]]$d$x[grids[[i]]$m])
          
          grids[[i]]$d$p_a[grids[[i]]$m] <- 0
          
          grids[[i + 1]]$d$p_a[1] <- 0
        } else {
          grids[[i + 1]]$d$x[1] <- grids[[i]]$d$x[grids[[i]]$m + 1]
          
          grids[[i]]$d <- head(grids[[i]]$d, -1)
          
          
          grids[[i + 1]]$d$s_upper[1] <- a / (grids[[i + 1]]$d$x[2] - grids[[i + 1]]$d$x[1])
          
          grids[[i + 1]]$d$p_a[1] <- 0
          
          grids[[i]]$d$p_a[grids[[i]]$m] <- 0
        }
      } else {
        grids[[i]]$m <- grids[[i]]$m + 1
        
        grids[[i]]$d$s_upper[grids[[i]]$m] <- a / (grids[[i + 1]]$d$x[1] - grids[[i]]$d$x[grids[[i]]$m])
        
        grids[[i]]$d$p_a[grids[[i]]$m] <- 0
      }
      
      m1 <- grids[[i]]$m
      
      m2 <- grids[[i + 1]]$m
      
      area[2] <- area[2] + m1 * a
      
      g_len[length(g_len) + 1] <- m1
      
      
      
      if (i == (mode_n - 1)) {
        final_grid <- rbind(final_grid, grids[[i]]$d, grids[[i + 1]]$d)
        
        area[2] <- area[2] + m2 * a
        
        g_len[length(g_len) + 1] <- m2
      } else {
        final_grid <- rbind(final_grid, grids[[i]]$d)
      }
    }
  } else {
    final_grid <- grids[[1]]$d
    
    area[2] <- grids[[1]]$m * a
    
    g_len[length(g_len) + 1] <- grids[[1]]$m
  }
  
  steps_number <- sum(g_len) # m
  
  x1 <- final_grid$x[1]
  
  xm <- final_grid$x[steps_number + 1]
  
  area[1] <- cdf(x1)
  
  area[3] <- 1 - cdf(xm)
  
  normalizing_con <- sum(area)
  
  area_cum_sum <- cumsum(area)
  
  sampling_probabilities <- (area_cum_sum / normalizing_con)[1:2]
  
  unif_scaler <- normalizing_con / area[2]
  
  lt_properties <- rep(0,5)
  rt_properties <- rep(0,6)
  
  invisible(list(grid_data = final_grid, areas = area, steps_number = steps_number, sampling_probabilities = sampling_probabilities, unif_scaler = unif_scaler, lt_properties = lt_properties, rt_properties = rt_properties))
}







#' bi Steps Builder
#'
#' @param lb scaler density lower bound
#' @param rb scaler density upper bound
#' @param a scaler step area
#' @param th scaler pre_acceptance threshold
#' @param mode vector density modes
#' @param mode_i mode index
#' @param mode_n total number of modes
#' @param f density function
#' @param h log transform of the density function
#' @param h_prime first derivative of h
#'
#' @return
#'
#' @examples
bi_find_steps <- function(lb = -Inf, rb = Inf, a, th, mode, mode_i, mode_n, f) {
  
  memory_res <- (max(500, ceiling(1 / a)) + 500) / 2
  
  x <- rep(NA, memory_res * 2 + 1)
  s_upper <- rep(NA, memory_res * 2 + 1)
  s_lower <- rep(NA, memory_res * 2 + 1)
  p_a <- rep(NA, memory_res * 2 + 1)
  s_upper_lower <- rep(NA, memory_res * 2 + 1)
  
  
  r <- 0
  l <- 0
  
  i <- memory_res + 1
  
  
  if (mode != rb) {
    x_c <- mode
    
    f_x <- f(x_c)
    
    while (TRUE) {
      x_next <- x_c + a / f_x
      
      if (x_next > rb || (mode_i == mode_n && f(x_next) / f_x < th)) {
        x[i + r] <- x_c
        break
      }
      
      f_x_next <- f(x_next)
      
      if (f_x_next > f_x) {
        x[i + r] <- x_c
        s_upper[i + r] <- f_x
        r_tail_area <- 0
        break
      }
      
      
      x[i + r] <- x_c
      s_upper[i + r] <- f_x
      s_lower[i + r] <- f_x_next
      s_upper_lower[i + r] <- s_upper[i + r] / s_lower[i + r]
      p_a[i + r] <- s_lower[i + r] / s_upper[i + r]
      
      
      f_x <- f_x_next
      x_c <- x_next
      
      r <- r + 1
    }
  }
  
  
  
  if (mode != lb) {
    x_previous <- mode
    f_x_previous <- f(x_previous)
    
    
    while (TRUE) {
      
      x_c <- x_previous - a / f_x_previous
      
      
      if (x_c < lb || (mode_i == 1 && f(x_c) / f_x_previous  < th)) {
        break
      }
      
      f_x <- f(x_c)
      
      if (f(x_c) > f_x_previous) {
        break
      }
      
      l <- l + 1
      x[i - l] <- x_c
      s_upper[i - l] <- f_x_previous
      s_lower[i - l] <- f_x
      s_upper_lower[i - l] <- s_upper[i - l] / s_lower[i - l]
      p_a[i - l] <- f_x / f_x_previous
    
      f_x_previous <- f_x
      x_previous <- x_c
    }
  }
  
  m <- l + r
  
  d <- data.frame(
    x = x, s_upper = s_upper, p_a = p_a,
    s_upper_lower = s_upper_lower
  )
  
  d <- subset(d, rowSums(is.na(d)) != ncol(d))
  
  return(list(d = d, m = m))
}



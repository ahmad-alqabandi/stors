
#' @method plot grid
#' @export
plot.grid <- function(grid, x_min = NA, x_max = NA,...){
  
  
  n = nrow(grid$grid_data)
  
  f <- eval(parse(text = grid$dens_func))
  
  if(grid$tails_method == "ARS"){
    lf  <- function(x) exp( grid$lt_properties[5] * (x - grid$grid_data$x[1]) + grid$lt_properties[3] )
    
    rf  <- function(x) exp( grid$rt_properties[5] * (x - grid$grid_data$x[n]) + grid$rt_properties[6] )
  } else {
    
    lf <- f
    rf <- f
  }
  
  if(is.finite(grid$grid_bounds[1])){
    l_limit = 0
  }else{
    l_limit = 1
  }
  
  if(is.finite(grid$grid_bounds[2])){
    r_limit = 0
  }else{
    r_limit = 1
  }
    
  
  xx <- seq(from = grid$grid_data$x[1]-l_limit, to = grid$grid_data$x[n]+r_limit, by = min(0.01,grid$alpha))
  
  xl <- seq( from = grid$grid_data$x[1]-l_limit , to= grid$grid_data$x[1], by = 0.01)
  
  xr <- seq( from = grid$grid_data$x[n] , to= grid$grid_data$x[n]+r_limit, by = 0.01)
  
  xs = c(xl, xx,xr)
  
  yy <- f(xx)
  
  yl <- lf(xl)
  
  yr <- rf(xr)
  
  ys = c(yl, yy, yr)
  

  if(is.na(x_max))  x_max = xr[length(xr)]
  if(is.na(x_min))  x_min = xl[1]
  
  y_max = max(ys[ x_min <= xs &  xs <= x_max ])
  y_min = min(ys[x_min <= xs &  xs <= x_max])
  
  grid$grid_data$s_upper[n] = yr[1]
  
  grid$grid_data = rbind( c(grid$grid_data$x[1], yl[length(yl)], NA, NA), grid$grid_data)
  
  n <- n+1
  
  if (requireNamespace("ggplot2")) {
    suppressWarnings({
    ggplot2::ggplot() +
      ggplot2::geom_step(ggplot2::aes(grid$grid_data[1:(n),]$x, grid$grid_data[1:(n),]$s_upper), color = "red") +
      ggplot2::geom_step(ggplot2::aes(grid$grid_data[2:(n-1),]$x, grid$grid_data[2:(n-1),]$s_upper * grid$grid_data[2:(n-1),]$p_a), color = "green") +
      ggplot2::geom_line(ggplot2::aes(x = xx, y = yy), color = "black")  +
      ggplot2::geom_line(ggplot2::aes(x = xl, y = yl), color = "red") +
      ggplot2::geom_line(ggplot2::aes(x = xr, y = yr), color = "red") +
      ggplot2::geom_segment(ggplot2::aes(x=grid$grid_data$x[(n-1)], y=grid$grid_data$s_upper[(n-1)] * grid$grid_data$p_a[(n-1)], xend = grid$grid_data$x[(n)], yend = grid$grid_data$s_upper[(n-1)] * grid$grid_data$p_a[(n-1)]), color="green")+
      ggplot2::xlab("x") +
      ggplot2::ylab("Density") +
      ggplot2::coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))
    })
  } else{
    class(grid) <- "list"
    plot(grid, ...)
  }
  
}




#' Print All Grids
#'
#' @description
#' Prints details of all grids stored by the user, including grid name, size, efficiency, etc.
#'
#' @export
#'
#' @examples
#' # To print details of all stored grids
#' print_grids()
print_grids <- function( ){
  
  stopifnot(" there are no grids stored by the user " = nrow(stors_env$grids$user) != 0)
  
  grids <- list.files(path = stors_env$user_dirs$data_dir, full.names = TRUE)
  
  grids_details <- file.info(grids)
  
  grids_sizes <- grids_details[grids_details$isdir == FALSE,]$size
  
  print(stors_env$grids$user)
  
  cat("\n grids_size : ", sum(as.double(grids_sizes))/1028," KB")
  
}




#' @method print grid
#' @export
print.grid <-function(grid, ...){
  
  cat("The grid has ", grid$steps_number ," steps, in the domain range [", grid$grid_data$x[1] ,",", grid$grid_data$x[grid$steps_number+1],"].\n")
  cat(sprintf("With a sampling efficiency of %.2f%%", 1/sum(grid$areas) * 100))
  
}


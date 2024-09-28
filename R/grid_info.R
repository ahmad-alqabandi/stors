
#' Plot Method for Grid Objects
#'
#' This method extends the generic `plot` function for objects of class `grid`.
#' It offers custom plotting functionality specifically designed for visualizing grid objects.
#'
#' @description
#' This function evaluates the properties of the included target and proposal functions to create a plot for both functions. In cases where the
#' proposal function's steps part is too dense, `x_min` and `x_max` can be set to crop and scale the chart for better visualization.
#'
#' @param x A list generated using STORS' `build_grid()` or `grid_optimizer()` functions.
#' @param x_min A scalar that represents the left cropping of the chart on the x-axis.
#' @param x_max A scalar that represents the right cropping of the chart on the x-axis.
#' @param ... Additional arguments passed to the `plot` function.
#' 
#' @return A plot of the target density and proposal. If `ggplot2` is available, it returns a `ggplot` object representing the plot. otherwise, it uses the base `plot` function.
#'
#' @seealso
#' \code{\link{print.grid}}
#' 
#' 
#' @examples
#' # Define the density function, its logarithm,
#' # and its derivative for the standard normal distribution
#' modes_norm = 0
#' f_norm <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x^2) }
#' h_norm <- function(x) { log(f_norm(x)) }
#' h_prime_norm <- function(x) { -x }
#'
#' # Build a dense grid for the standard normal distribution
#' norm_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_norm,
#'  f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 4000)
#'
#' # Plot the generated grid
#' plot(norm_grid)
#'
#' # To visualize the grid in a cropped area between -0.1 and 0
#' plot(norm_grid, x_min = -0.1, x_max = 0)
#' 
#'@method plot grid
#' @export
plot.grid <- function(x, x_min = NA, x_max = NA,...){
  
  grid <- x
  
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
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    #suppressWarnings({
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
    #})
  } else{
    class(grid) <- "list"
    plot(grid, ...)
  }
  
}




#' Print Method for Grid Objects
#'
#' This method extends the generic `print` function for objects of class `grid`.
#' It prints the provided grid's features such as the number of steps, steps limit, and efficiency.
#'
#' @description
#' The function displays detailed information about the grid object created by STORS' `build_grid()` or `grid_optimizer()` functions.
#'  This includes the number of steps within the grid, the range of values covered by the grid, and the grid's sampling efficiency.
#'   This information is crucial for understanding the structure and performance of the grid in sampling processes.
#'
#' @param x A list generated using STORS' `build_grid()` or `grid_optimizer()` functions.
#' @param ... Additional arguments passed to the `print` function.
#' 
#' @return Prints a summary of the grid's properties, but does not return any value.
#'
#' @examples
#' # Define the density function, its logarithm,
#' #and its derivative for the standard normal distribution
#' modes_norm = 0
#' f_norm <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x^2) }
#' h_norm <- function(x) { log(f_norm(x)) }
#' h_prime_norm <- function(x) { -x }
#'
#' # Build a dense grid for the standard normal distribution
#' norm_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_norm,
#'  f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)
#'
#' # Print the properties of the generated grid
#' print(norm_grid)
#'
#' @export
#' @method print grid
print.grid <- function(x, ...) {
  
  grid <- x
  
  formatted_steps <- format(grid$steps_number, big.mark = ",", scientific = FALSE)
  cat("The grid contains", formatted_steps, "steps within the domain range  [", grid$grid_data$x[1], ",", grid$grid_data$x[grid$steps_number + 1], "].\n")
  cat(sprintf("With a sampling efficiency of %.2f%%", 1 / sum(grid$areas) * 100), "\n")
}




#' Print Grids
#'
#' @description
#' This function prints details of all grids stored by the user. It provides information on each grid, including the grid name, size, efficiency, and other relevant details.
#'
#' @export
#'
#' @examples
#' # First, let's create a grid to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_grid = build_grid(f = f_normal, modes = 0, lb = -Inf, rb = Inf, steps = 1000)
#' print(normal_grid)
#'
#' # `print_grids()` prints all grids stored in R's internal data directory.
#' # To see this, we first save 'normal_grid' using `save_grid()`
#' save_grid(normal_grid, "normal")
#'
#' # Since 'normal_grid' is now stored on this machine,
#' # we can confirm this by printing all saved grids
#' print_grids()
#' 
print_grids <- function( ){
  
  if(nrow(stors_env$grids$user) == 0){
    message("No grids are currently stored.")
  }else{
    
    grids <- list.files(path = stors_env$user_dirs$data_dir, full.names = TRUE)
    
    grids_details <- file.info(grids)
    
    grids_sizes <- grids_details[grids_details$isdir == FALSE,]$size
    
    print(stors_env$grids$user)
    
    grids_sizes <- formatC(sum(as.double(grids_sizes))/1028, format = "f", digits = 2)
    
    cat("\n grids_size : ", grids_sizes," KB")
  }

}







#' Save User Grid
#'
#' @description
#' This function stores grids generated by the `build_grid()` function in R's internal data directory. It is useful when users want to reuse a grid across multiple R sessions.
#'
#' @param grid list representing an optimized grid generated using the `build_grid()` function.
#' @param grid_name string specifying the name under which the grid will be saved.
#'
#' @return
#' This function will produce an error if the grid is not generated by the `build_grid()` function. Otherwise, it successfully saves the grid without returning any value upon completion.
#' @export
#' 
#' 
#' @import digest digest
#' 
#' @examples
#' # First, let's create a grid to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_grid = build_grid(f = f_normal, modes = 0, lb = -Inf, rb = Inf, steps = 1000)
#' print(normal_grid)
#'
#' # Then, we can save this grid in R's internal data directory using `save_grid()`
#' # with the name "normal"
#' save_grid(normal_grid, "normal")
#'
#' # To make sure the 'normal_grid' has been stored in R's internal data directory,
#' # we can print all saved grids using `print_grids()`
#' print_grids()
#' 
save_grid <- function(grid, grid_name) {
  
  is_valid_grid(grid)
  
  if(grid_name %in% stors_env$grids$user$name) invisible(delete_grid(grid_name))
  
  grids_file_path <- file.path(stors_env$user_dirs$data_dir, paste0(grid_name, ".rds"))
  
  saveRDS(grid, grids_file_path)
  
  efficiency <- (1/sum(grid$areas))
  
  stors_env$grids$user[ nrow(stors_env$grids$user) + 1 ,] = list( grid_name, efficiency)
}


#' Delete Grid
#'
#' @description
#' This function deletes a grid that was previously stored by the user using the `save_grid()` function. It is useful for managing stored grids and freeing up space.
#'
#' @param grid_name A string specifying the name of the grid to be deleted.
#'
#' @return
#' If `grid_name` does not exist, the function returns an error message. If the grid exists and is successfully deleted, a message confirming its successful removal will be displayed.
#'
#' @export
#'
#' @examples
#' # First, let's create a grid to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_grid = build_grid(f = f_normal, modes = 0, lb = -Inf, rb = Inf, steps = 1000)
#' print(normal_grid)
#'
#' # Then, save this grid in R's internal data directory using `save_grid()` with the name "normal"
#' save_grid(normal_grid, "normal")
#'
#' # Now, we can print all grids stored on this machine using `print_grids()`
#' print_grids()
#'
#' # The list will include the 'normal_grid' stored under the name "normal"
#'
#' # To delete the "normal" grid from the machine, pass its name to `delete_grid`
#' delete_grid("normal")
#'
#' # Now, when we print all stored grids, the "normal" grid will no longer be listed
#' print_grids()
#' 
delete_grid <- function(grid_name){
  
  stopifnot("This grid does not exist." = grid_name %in% stors_env$grids$user$name)
  
  file.remove(file.path( stors_env$user_dirs$data_dir, paste0(grid_name,".rds") ))
  stors_env$grids$user = stors_env$grids$user[stors_env$grids$user$name != grid_name,]
  cat(grid_name, "grid has been deleted successfully")
  
}


#' Load Stored Grid
#'
#' @description
#' This function loads a grid into memory that was previously saved using the `save_grid()` function. It is useful for retrieving saved grids for further analysis or processing.
#'
#' @param grid_name A string specifying the name of the grid to be loaded.
#'
#' @return
#' Returns a list representing the grid stored under `grid_name` in R's internal data directory. If the grid corresponding to the specified name does not exist, an error message is displayed.
#'
#' @export
#'
#' @examples
#' # First, let's create a grid to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_grid = build_grid(f = f_normal, modes = 0, lb = -Inf, rb = Inf, steps = 1000)
#' print(normal_grid)
#'
#' # Then, save this grid in R's internal data directory using `save_grid()` with the name "normal"
#' save_grid(normal_grid, "normal")
#'
#' # Now, in case the R session is restarted and the grid is no longer in memory,
#' # it can be loaded from the machine as follows:
#' loaded_normal_grid <- load_grid("normal")
#' print(loaded_normal_grid)
#'  
load_grid <- function(grid_name) {
  
  if (!(grid_name %in% stors_env$grids$user$name)) {
    stop("There is no grid named '", grid_name, "' stored on your machine.")
  }
  
  grids_file_path <- file.path(stors_env$user_dirs$data_dir, paste0(grid_name, ".rds"))
  
  grid <- readRDS(grids_file_path)
  
  stors_env$created_girds_Id  = append(stors_env$created_girds_Id , digest(grid))
  
  return(grid)
  
}




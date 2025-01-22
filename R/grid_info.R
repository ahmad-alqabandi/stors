
#' Plot Method for Grid Objects
#'
#' This method extends the generic \code{plot()} function for objects of class \code{grid}.
#' It offers custom plotting functionality specifically designed for visualizing grid objects.
#'
#' @description
#' This function evaluates the properties of the included target and proposal functions to create a plot for both functions. In cases where the
#' proposal function's steps part is too dense, \code{x_min} and \code{x_max} can be set to crop and scale the chart for better visualization.
#'
#' @param x A list generated using STORS' \code{build_grid()} or \code{grid_optimizer()} functions.
#' @param x_min A scalar that represents the left cropping of the chart on the x-axis.
#' @param x_max A scalar that represents the right cropping of the chart on the x-axis.
#' @param ... Additional arguments passed to the \code{plot} function.
#' @return A plot of the target density and proposal. If \code{ggplot2} is available, it returns a \code{ggplot} object representing the plot. otherwise, it uses the base \code{plot()} function.
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
#' @method plot grid
#' @export
plot.grid <- function(x,
                      x_min = NA,
                      x_max = NA,
                      ...) {
  grid <- x
  n <- nrow(grid$grid_data)
  f <- eval(parse(text = grid$dens_func))
  f <- create_function(f, grid$density_arguments)
  lf <- f
  rf <- f



  if (grid$tails_method == "ARS") {
    if (!all(grid$lt_properties == 0)) {
      lf  <- function(x) exp(grid$lt_properties[5] * (x - grid$grid_data$x[1]) + grid$lt_properties[3])
    }

    if (!all(grid$rt_properties == 0)) {
      rf  <- function(x) exp(grid$rt_properties[5] * (x - grid$grid_data$x[n]) + grid$rt_properties[6])
    }
  }


  if (is.finite(grid$grid_bounds[1])) {
    x_from <- grid$grid_bounds[1]
  } else {
    x_from <- grid$grid_data$x[1] - 5
  }

  if (is.finite(grid$grid_bounds[2])) {
    x_to <- grid$grid_bounds[2]
  } else {
   x_to <- grid$grid_data$x[n] + 5
  }

  xx <- seq(
    from = x_from,
    to = x_to,
    by = min(0.01, grid$alpha)
  )
  xl <- seq(
    from = x_from,
    to = grid$grid_data$x[1],
    by = 0.01
  )
  xr <- seq(
    from = grid$grid_data$x[n],
    to = x_to,
    by = 0.01
  )
  xs <- c(xl, xx, xr)

  yy <- f(xx)
  yl <- lf(xl)
  yr <- rf(xr)
  ys <- c(yl, yy, yr)

  if (is.na(x_max))
    x_max <- xr[length(xr)]

  if (is.na(x_min))
    x_min <- xl[1]

  y_max <- max(ys[x_min <= xs &  xs <= x_max])
  y_min <- min(ys[x_min <= xs &  xs <= x_max])

  grid$grid_data$s_upper[n] <- yr[1]
  grid$grid_data <- rbind(c(grid$grid_data$x[1], yl[length(yl)], NA, NA), grid$grid_data)

  n <- n + 1

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressWarnings({
    ggplot2::ggplot() +
      ggplot2::geom_step(ggplot2::aes(grid$grid_data[1:(n), ]$x, grid$grid_data[1:(n), ]$s_upper),
                         color = "red") +
      ggplot2::geom_step(
        ggplot2::aes(
          grid$grid_data[2:(n - 1), ]$x,
          grid$grid_data[2:(n - 1), ]$s_upper * grid$grid_data[2:(n - 1), ]$p_a
        ),
        color = "green"
      ) +
      ggplot2::geom_line(ggplot2::aes(x = xx, y = yy), color = "black")  +
      ggplot2::geom_line(ggplot2::aes(x = xl, y = yl), color = "red", linetype = "dotted") +
      ggplot2::geom_line(ggplot2::aes(x = xr, y = yr), color = "red", linetype = "dotted") +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = grid$grid_data$x[(n - 1)],
          y = grid$grid_data$s_upper[(n - 1)] * grid$grid_data$p_a[(n - 1)],
          xend = grid$grid_data$x[(n)],
          yend = grid$grid_data$s_upper[(n - 1)] * grid$grid_data$p_a[(n - 1)]
        ),
        color = "green"
      ) +
      ggplot2::xlab("x") +
      ggplot2::ylab("Density") +
      ggplot2::coord_cartesian(xlim = c(x_min, x_max),
                               ylim = c(y_min, y_max))
    })
  } else {
    class(grid) <- "list"
    plot(grid, ...)
  }

}


#' Print Method for Grid Objects
#'
#' This method extends the generic \code{print} function for objects of class \code{grid}.
#' It prints the provided grid's features such as the number of steps, steps limit, and efficiency.
#'
#' @description
#' The function displays detailed information about the grid object created by STORS' \code{build_grid()} or \code{grid_optimizer()} functions.
#'  This includes the number of steps within the grid, the range of values covered by the grid, and the grid's sampling efficiency.
#'   This information is crucial for understanding the structure and performance of the grid in sampling processes.
#'
#' @param x A list generated using STORS' \code{build_grid()} or \code{grid_optimizer()} functions.
#' @param ... Additional arguments passed to the \code{print} function.
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
#'
#' print(norm_grid)
#' @method print grid
#' @export
print.grid <- function(x, ...) {

  grid <- x

  # Format the number of steps with commas and without scientific notation
  formatted_steps <- format(grid$steps_number, big.mark = ",", scientific = FALSE)

  # Calculate the domain range
  domain_start <- grid$grid_data$x[1]
  n <- length(grid$grid_data$x)
  domain_end <- grid$grid_data$x[n]

  # Calculate sampling efficiency
  sampling_efficiency <- (grid$target_function_area / (sum(grid$areas))) * 100

  # Improved printing with clearer structure
  message("\n=========================================\n")
  message("Grid Summary\n")
  message("-----------------------------------------\n")
  message("Total steps:      ", formatted_steps, "\n")
  message("Steps range:     [", sprintf("%.6f", domain_start), ", ", sprintf("%.6f", domain_end), "]\n")
  message(sprintf("Sampling efficiency: %.2f%%", sampling_efficiency), "\n")
  message("=========================================\n\n")

}


#' Print Grids
#'
#' @description
#' This function prints details of all grids stored by the user. It provides information on each grid, including the grid name, size, efficiency, and other relevant details.
#'
#'
#'
#' @examples
#' # First, let's create a grid to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_grid = build_grid(f = f_normal, modes = 0, lb = -Inf, rb = Inf, steps = 1000)
#' print(normal_grid)
#'
#' # `print_grids()` prints all grids stored in R's internal data directory.
#' # To see this, we first save `normal_grid` using `save_grid()`
#' save_grid(normal_grid, "normal")
#'
#' # Since `normal_grid` is now stored on this machine,
#' # we can confirm this by printing all saved grids
#' print_grids()
#'
#' # Example 2: Create and Save a Grid for a Bimodal Distribution
#' f_bimodal <- function(x) {
#'   0.5 * (1 / sqrt(2 * pi)) * exp(-(x^2) / 2) +
#'   0.5 * (1 / sqrt(2 * pi)) * exp(-((x - 4)^2) / 2)
#' }
#' modes_bimodal = c(0, 4)
#' bimodal_grid = build_grid(f = f_bimodal, modes = modes_bimodal, lb = -Inf, rb = Inf, steps = 1000)
#' save_grid(bimodal_grid, "bimodal")
#' print(bimodal_grid)
#'
#' # To print all stored grids after saving bimodal_grid
#' print_grids()
#'
#' @export
print_grids <- function() {

  user_grids <- list.files(stors_env$user_grids_dir)

  if (length(user_grids) == 0) {
    message("No grids are currently stored.")
  } else {
    grids <- list.files(path = stors_env$user_grids_dir,
                        full.names = TRUE)
    grids_details <- file.info(grids)
    grids_sizes <- grids_details[grids_details$isdir == FALSE, ]$size
    grids_sizes <- formatC(sum(as.double(grids_sizes)) / 1028,
                           format = "f",
                           digits = 2)
    message("grids_size : ", grids_sizes, " KB")
  }

}







#' Save User Grid
#'
#' @description
#' This function stores grids generated by the \code{build_grid()} function in R's internal data directory. It is useful when users want to reuse a grid across multiple R sessions.
#'
#' @param grid list representing an optimized grid generated using the \code{build_grid()} function.
#' @param grid_name string specifying the name under which the grid will be saved.
#'
#' @return
#' This function will produce an error if the grid is not generated by the \code{build_grid()} function. Otherwise, it successfully saves the grid without returning any value upon completion.
#'
#' @examples
#' # First, let's create a grid to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp( -0.5 * x^2) }
#' normal_grid = build_grid(f = f_normal, modes = 0, lb = -Inf, rb = Inf, steps = 1000)
#' print(normal_grid)
#'
#' # Then, we can save this grid in R's internal data directory using `save_grid()`
#' # with the name "normal"
#' save_grid(normal_grid, "normal")
#'
#' # To make sure the `normal_grid` has been stored in R's internal data directory,
#' # we can print all saved grids using `print_grids()`
#' print_grids()
#'
#'
#' @import digest digest
#' @export
save_grid <- function(grid, grid_name) {
  if (!is_valid_grid(grid))
    stop("This grid is not valid")

  grids_file_path <- file.path(stors_env$user_grids_dir, paste0(grid_name, ".rds"))
  saveRDS(grid, grids_file_path)
}


#' Delete Grid
#'
#' @description
#' This function deletes a grid that was previously stored by the user using the \code{save_grid()} function. It is useful for managing stored grids and freeing up space.
#'
#' @param grid_name A string specifying the name of the grid to be deleted.
#'
#' @return
#' If \code{grid_name} does not exist, the function returns an error message. If the grid exists and is successfully deleted, a message confirming its successful removal will be displayed.
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
#' # The list will include the `normal_grid` stored under the name "normal"
#'
#' # To delete the "normal" grid from the machine, pass its name to `delete_grid`
#' delete_grid("normal")
#'
#' # Now, when we print all stored grids, the "normal" grid will no longer be listed
#' print_grids()
#'
delete_grid <- function(grid_name) {

  user_grids <- list.files(stors_env$user_grids_dir)
  grid_name <- paste0(grid_name, ".rds")

  stopifnot("This grid does not exist." = grid_name %in% user_grids)

  file.remove(file.path(stors_env$user_grids_dir, grid_name))
  message(grid_name, "grid deleted successfully")

}


#' Load Stored Grid
#'
#' @description
#' This function loads a grid into memory that was previously saved using the \code{save_grid()} function. It is useful for retrieving saved grids for further analysis or processing.
#'
#' @param grid_name A string specifying the name of the grid to be loaded.
#'
#' @return
#' Returns a list representing the grid stored under \code{grid_name} in R's internal data directory. If the grid corresponding to the specified name does not exist, an error message is displayed.
#'
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
#' @export
load_grid <- function(grid_name) {

  grid_name <- paste0(grid_name, ".rds")

  user_grids <- list.files(stors_env$user_grids_dir)

  if (grid_name %in% user_grids) {

    grid_path <- file.path(stors_env$user_grids_dir, grid_name)

    grid <- readRDS(grid_path)

    if (!is_valid_grid(grid))
      stop("This grid is not valid")

    return(grid)


  }else {
    stop("There is no grid named '",
         grid_name,
         "' stored on your machine.")
  }


}

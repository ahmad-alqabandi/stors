#' Plot Method for Proposal Objects
#'
#' This method extends the generic \code{plot()} function for objects of class \code{proposal}.
#' It offers custom plotting functionality specifically designed for visualizing proposal objects.
#'
#' @description
#' This function evaluates the properties of the included target and proposal functions to create a plot for both functions. In cases where the
#' proposal function's steps part is too dense, \code{x_min} and \code{x_max} can be set to crop and scale the chart for better visualization.
#'
#' @param x A list generated using STORS' \code{build_proposal()} or \code{proposal_optimizer()} functions.
#' @param x_min A scalar that represents the left cropping of the chart on the x-axis.
#' @param x_max A scalar that represents the right cropping of the chart on the x-axis.
#' @param ... Additional arguments passed to the \code{plot} function.
#' @return A plot of the target density and proposal. If \code{ggplot2} is available, it returns a \code{ggplot} object representing the plot. otherwise, it uses the base \code{plot()} function.
#'
#' @seealso
#' \code{\link{print.proposal}}
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
#' # Build a dense proposal for the standard normal distribution
#' norm_proposal = build_proposal(lower = -Inf, upper = Inf, mode = modes_norm,
#'  f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 4000)
#'
#' # Plot the generated proposal
#' plot(norm_proposal)
#'
#' # To visualize the proposal in a cropped area between -0.1 and 0
#' plot(norm_proposal, x_min = -0.1, x_max = 0)
#'
#' @method plot proposal
#' @export
plot.proposal <- function(x,
                          x_min = NA,
                          x_max = NA,
                          ...) {
  proposal <- x
  n <- nrow(proposal$data)
  f <- eval(parse(text = proposal$dens_func))
  f <- create_function(f, proposal$density_arguments)
  lf <- f
  rf <- f



  if (proposal$tails_method == "ARS") {
    if (!all(proposal$lt_properties == 0)) {
      lf  <- function(x) exp(proposal$lt_properties[5] * (x - proposal$data$x[1]) + proposal$lt_properties[3])
    }

    if (!all(proposal$rt_properties == 0)) {
      rf  <- function(x) exp(proposal$rt_properties[5] * (x - proposal$data$x[n]) + proposal$rt_properties[6])
    }
  }


  if (is.finite(proposal$proposal_bounds[1])) {
    x_from <- proposal$proposal_bounds[1]
  } else {
    x_from <- proposal$data$x[1] - 5
  }

  if (is.finite(proposal$proposal_bounds[2])) {
    x_to <- proposal$proposal_bounds[2]
  } else {
    x_to <- proposal$data$x[n] + 5
  }

  xx <- seq(from = x_from,
            to = x_to,
            by = min(0.01, proposal$alpha))
  xl <- seq(from = x_from,
            to = proposal$data$x[1],
            by = 0.01)
  xr <- seq(from = proposal$data$x[n],
            to = x_to,
            by = 0.01)
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

  proposal$data$s_upper[n] <- yr[1]
  proposal$data <- rbind(c(proposal$data$x[1], yl[length(yl)], NA, NA), proposal$data)

  n <- n + 1

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressWarnings({
      ggplot2::ggplot() +
        ggplot2::geom_step(ggplot2::aes(proposal$data[1:(n), ]$x, proposal$data[1:(n), ]$s_upper),
                           color = "red") +
        ggplot2::geom_step(
          ggplot2::aes(
            proposal$data[2:(n - 1), ]$x,
            proposal$data[2:(n - 1), ]$s_upper * proposal$data[2:(n - 1), ]$p_a
          ),
          color = "green"
        ) +
        ggplot2::geom_line(ggplot2::aes(x = xx, y = yy), color = "black")  +
        ggplot2::geom_line(ggplot2::aes(x = xl, y = yl),
                           color = "red",
                           linetype = "dotted") +
        ggplot2::geom_line(ggplot2::aes(x = xr, y = yr),
                           color = "red",
                           linetype = "dotted") +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = proposal$data$x[(n - 1)],
            y = proposal$data$s_upper[(n - 1)] * proposal$data$p_a[(n - 1)],
            xend = proposal$data$x[(n)],
            yend = proposal$data$s_upper[(n - 1)] * proposal$data$p_a[(n - 1)]
          ),
          color = "green"
        ) +
        ggplot2::xlab("x") +
        ggplot2::ylab("Density") +
        ggplot2::coord_cartesian(xlim = c(x_min, x_max),
                                 ylim = c(y_min, y_max))
    })
  } else {
    class(proposal) <- "list"
    plot(proposal, ...)
  }

}


#' Print Method for proposal Objects
#'
#' This method extends the generic \code{print} function for objects of class \code{proposal}.
#' It prints the provided proposal's features such as the number of steps, steps limit, and efficiency.
#'
#' @description
#' The function displays detailed information about the proposal object created by STORS' \code{build_proposal()} or \code{proposal_optimizer()} functions.
#'  This includes the number of steps within the proposal, the range of values covered by the proposal, and the proposal's sampling efficiency.
#'   This information is crucial for understanding the structure and performance of the proposal in sampling processes.
#'
#' @param x A list generated using STORS' \code{build_proposal()} or \code{proposal_optimizer()} functions.
#' @param ... Additional arguments passed to the \code{print} function.
#'
#' @return Prints a summary of the proposal's properties, but does not return any value.
#'
#' @examples
#' # Define the density function, its logarithm,
#' #and its derivative for the standard normal distribution
#' modes_norm = 0
#' f_norm <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x^2) }
#' h_norm <- function(x) { log(f_norm(x)) }
#' h_prime_norm <- function(x) { -x }
#'
#' # Build a dense proposal for the standard normal distribution
#' norm_proposal = build_proposal(lower = -Inf, upper = Inf, mode = modes_norm,
#'  f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)
#'
#' # Print the properties of the generated proposal
#'
#' print(norm_proposal)
#' @method print proposal
#' @export
print.proposal <- function(x, ...) {
  proposal <- x

  formatted_steps <- format(proposal$steps_number,
                            big.mark = ",",
                            scientific = FALSE)

  # Calculate the domain range
  domain_start <- proposal$data$x[1]
  n <- length(proposal$data$x)
  domain_end <- proposal$data$x[n]

  # Calculate sampling efficiency
  sampling_efficiency <- (proposal$target_function_area / sum(proposal$areas)) * 100

  # Clearer formatting using `cli`
  cli::cli_h1("Proposal Summary")

  cli::cli_ul(c(
    "Total steps: {formatted_steps}",
    "Steps range: [{sprintf('%.6f', domain_start)}, {sprintf('%.6f', domain_end)}]",
    "Sampling efficiency: {sprintf('%.2f%%', sampling_efficiency)}"
  ))
}


#' Print proposals
#'
#' @description
#' This function prints details of all proposals stored by the user. It provides information on each proposal, including the proposal name, size, efficiency, and other relevant details.
#'
#' @return Prints a summary of all proposals, but does not return any value.
#'
#' @examples
#' # First, let's create a proposal to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_proposal = build_proposal(f = f_normal, modes = 0, lower = -Inf, upper = Inf, steps = 1000)
#' print(normal_proposal)
#'
#' # `print_proposals()` prints all proposals stored in R's internal data directory.
#' # To see this, we first save `normal_proposal` using `save_proposal()`
#' save_proposal(normal_proposal, "normal")
#'
#' # Since `normal_proposal` is now stored on this machine,
#' # we can confirm this by printing all saved proposals
#' print_proposals()
#'
#' # Example 2: Create and Save a proposal for a Bimodal Distribution
#' f_bimodal <- function(x) {
#'   0.5 * (1 / sqrt(2 * pi)) * exp(-(x^2) / 2) +
#'   0.5 * (1 / sqrt(2 * pi)) * exp(-((x - 4)^2) / 2)
#' }
#' modes_bimodal = c(0, 4)
#' bimodal_proposal = build_proposal(f = f_bimodal, modes = modes_bimodal,
#' lower = -Inf, upper = Inf, steps = 1000)
#'
#' save_proposal(bimodal_proposal, "bimodal")
#' print(bimodal_proposal)
#'
#' # To print all stored proposals after saving bimodal_proposal
#' print_proposals()
#'
#' @export
print_proposals <- function() {
  user_proposals <- list.files(stors_env$user_proposals_dir)

  if (length(user_proposals) == 0) {
    cli::cli_inform("{.strong No stored proposals.}")
  } else {
    files <- list.files(path = stors_env$user_proposals_dir, full.names = TRUE)
    files_details <- file.info(files)

    files_names <- tools::file_path_sans_ext(basename(files))
    files_date <- file.mtime(files)
    files_sizes <- files_details[files_details$isdir == FALSE, ]$size

    files_sizes_kb <- formatC(as.double(files_sizes) / 1028, format = "f", digits = 2)
    files_total_size_kb <- formatC(as.double(sum(files_sizes)) / 1028, format = "f", digits = 2)

    cli::cli_h1("Proposals Data")

    name_width <- 20
    size_width <- 15
    date_width <- 20

    cli::cli_verbatim(sprintf("%-*s | %-*s | %-*s",
                              name_width, "Name",
                              size_width, "Size (KB)",
                              date_width, "Date"))
    cat("-----------------------------------------------\n")

    for (k in seq_along(files_names)) {
      cli::cli_verbatim(sprintf("%-*s | %-*s | %-*s",
                                name_width, files_names[k],
                                size_width, paste0(files_sizes_kb[k]),
                                date_width, format(files_date[k], "%Y-%m-%d %H:%M:%S")))
    }

    cat("-----------------------------------------------\n")

    cli::cli_verbatim(sprintf("Total Size: %s", files_total_size_kb))

  }

}





#' Save User Proposal
#'
#' @description
#' This function stores proposals generated by the \code{build_proposal()} function in R's internal data directory. It is useful when users want to reuse a proposal across multiple R sessions.
#'
#' @param proposal list representing an optimized proposal generated using the \code{build_proposal()} function.
#' @param proposal_name string specifying the name under which the proposal will be saved.
#'
#' @return
#' This function will produce an error if the proposal is not generated by the \code{build_proposal()} function. Otherwise, it successfully saves the proposal without returning any value upon completion.
#'
#' @examples
#' # First, let's create a proposal to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp( -0.5 * x^2) }
#' normal_proposal = build_proposal(f = f_normal, modes = 0, lower = -Inf, upper = Inf, steps = 1000)
#' print(normal_proposal)
#'
#' # Then, we can save this proposal in R's internal data directory using `save_proposal()`
#' # with the name "normal"
#' save_proposal(normal_proposal, "normal")
#'
#' # To make sure the `normal_proposal` has been stored in R's internal data directory,
#' # we can print all saved proposals using `print_proposals()`
#' print_proposals()
#'
#'
#' @import digest digest
#' @export
save_proposal <- function(proposal, proposal_name) {

  if (!is_valid_proposal(proposal)) {
    cli::cli_abort(c("x" = "{.strong This proposal is not valid.}",
                     "i" = "Only proposals created using {.fn build_proposal} can be used."))
  }

  proposals_file_path <- file.path(stors_env$user_proposals_dir,
                                   paste0(proposal_name, ".rds"))
  saveRDS(proposal, proposals_file_path)
}


#' Delete Proposal
#'
#' @description
#' This function deletes a proposal that was previously stored by the user using the \code{save_proposal()} function. It is useful for managing stored proposals and freeing up space.
#'
#' @param proposal_name A string specifying the name of the proposal to be deleted.
#'
#' @return
#' If \code{proposal_name} does not exist, the function returns an error message. If the proposal exists and is successfully deleted, a message confirming its successful removal will be displayed.
#'
#' @export
#'
#' @examples
#' # First, let's create a proposal to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_proposal = build_proposal(f = f_normal, modes = 0, lower = -Inf, upper = Inf, steps = 1000)
#' print(normal_proposal)
#'
#' # Then, save this proposal in R's internal data directory using
#' # `save_proposal()` with the name "normal"
#' save_proposal(normal_proposal, "normal")
#'
#' # Now, we can print all proposals stored on this machine using `print_proposals()`
#' print_proposals()
#'
#' # The list will include the `normal_proposal` stored under the name "normal"
#'
#' # To delete the "normal" proposal from the machine, pass its name to `delete_proposal`
#' delete_proposal("normal")
#'
#' # Now, when we print all stored proposals, the "normal" proposal will no longer be listed
#' print_proposals()
#'
delete_proposal <- function(proposal_name) {
  user_proposals <- list.files(stors_env$user_proposals_dir)
  proposal_name <- paste0(proposal_name, ".rds")

  if (!(proposal_name %in% user_proposals)) {
    cli::cli_abort(c("x" = "The proposal {proposal_name} does not exist.",
                     "i" = "Please provide a valid proposal from the available proposals: {paste(user_proposals, collapse = ', ')}."))
  }

  file.remove(file.path(stors_env$user_proposals_dir, proposal_name))
  cli::cli_alert_success("{.val {proposal_name}} proposal deleted successfully")

}


#' Load Stored Proposal
#'
#' @description
#' This function loads a proposal into memory that was previously saved using the \code{save_proposal()} function. It is useful for retrieving saved proposals for further analysis or processing.
#'
#' @param proposal_name A string specifying the name of the proposal to be loaded.
#'
#' @return
#' Returns a list representing the proposal stored under \code{proposal_name} in R's internal data directory. If the proposal corresponding to the specified name does not exist, an error message is displayed.
#'
#' @examples
#' # First, let's create a proposal to sample from a standard normal distribution
#' f_normal <- function(x) { 0.3989423 * exp(-0.5 * x^2) }
#' normal_proposal = build_proposal(f = f_normal, modes = 0, lower = -Inf, upper = Inf, steps = 1000)
#' print(normal_proposal)
#'
#' # Then, save this proposal in R's internal data directory using
#' # `save_proposal()` with the name "normal"
#' save_proposal(normal_proposal, "normal")
#'
#' # Now, in case the R session is restarted and the proposal is no longer in memory,
#' # it can be loaded from the machine as follows:
#' loaded_normal_proposal <- load_proposal("normal")
#' print(loaded_normal_proposal)
#' @export
load_proposal <- function(proposal_name) {
  proposal_name <- paste0(proposal_name, ".rds")

  user_proposals <- list.files(stors_env$user_proposals_dir)

  if (proposal_name %in% user_proposals) {
    proposal_path <- file.path(stors_env$user_proposals_dir, proposal_name)

    proposal <- readRDS(proposal_path)

    if (!is_valid_proposal(proposal))
      cli::cli_abort(c("x" = "This proposal is not valid",
                       "i" = "Only proposals created using {.fn build_proposal} can be used."))

    return(proposal)

  } else {
    cli::cli_abort(c("x" = "There is no proposal named {.val {proposal_name}} stored on your machine",
                     "i" = "Expected location: {.path {proposal_path}}",
                     "i" = "Please check the name or ensure the file exists."))
  }
}

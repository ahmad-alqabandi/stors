#' Build User Grid
#'
#' @description
#' This function is essential for sampling from any uni-modal or multi-modal log-concave density function.
#'  It generates a proposal grid that can be utilized for this purpose.
#'  Simply provide the density function \code{f} and its modes. For enhanced accuracy, include
#' the log-transform \code{h} of the density and its first derivative \code{h_prime}.
#'  The grid optimization is pre-configured for efficiency,
#' but you can customize the grid-building process through various parameters to suit your specific needs.
#'
#' @param lb Scalar representing the lower bound of the target density.
#' @param rb Scalar representing the upper bound of the target density.
#' @param modes Vector indicating the modes of the density function.
#' @param f Function accepting a single argument, returning the probability density of the target.
#' @param h Function accepting a single argument, returning the log-transform of the target density.
#' @param h_prime Function accepting a single argument, returning the first derivative of the log-transformed target density.
#' @param steps Optional Scalar integer indicating the number of steps in the proposal distribution.
#' @param grid_range Optional Vector of two elements specifying the start and end points for constructing steps along the x-axis.
#' @param theta Optional Scalar defining the pre-acceptance threshold,
#'  dictating when the proposal steps constructing break based on the probability of pre-acceptance.
#' @param target_sample_size Scalar integer indicating the target sample size. The grid optimization process will take this number into account.
#' @param verbose Boolean if set to True, a table detailing the optimization areas and steps will be displayed during grid optimization.
#'
#' @details
#' The grid building process is executed through the construction of constant area rectangles,
#'  starting from the modes of the target distribution.
#' For each mode, rectangles are formed as steps around it,
#'  with a width defined by \eqn{(x_i - x_{i-1})} and a height determined by \eqn{\max(f(x_{i-1}), f(x_i))}.
#' This method effectively covers the target distribution in a stepped pattern.
#'
#'
#' The function \code{build_final_grid()} manages the construction of these steps and calculates values critical for the sampling process.
#'  When the resultant grid is used with the \code{stors()} function, these values are cached,
#' significantly enhancing the computational efficiency and hence improving sampling speed.
#'  During the optimization process, we aim for a certain grid
#' size based on L1-3 memory cache size. Therefore, we test the speed of grids of sizes \eqn{2^m} Kb.
#'  To achieve this, we estimate the uniform step area
#' based on a certain steps number that leads to the target cache size,
#'  \eqn{ \alpha = \frac{1}{\text{number of steps}} }.
#'
#'
#' The speed testing for each possible grid is initially based on a sample size of 1000.
#'  However, if the user wishes to optimize the grid for a different sample size, they can do so
#' by specifying the desired sample size using the \code{target_sample_size} argument.
#'
#'
#' In case the user wants to select a specific number of steps for the proposal grid
#' and bypass the optimization process, this can be done by specifying a steps number greater than the number of modes by 2 using the \code{steps} argument.
#'  If the target density is heavy-tailed,
#'   and the user wishes to stop the grid building process at a certain pre-acceptance threshold, this can be achieved by setting
#' the acceptance probability threshold \code{theta} \eqn{\theta}.
#'  Once the steps reach this level of pre-acceptance probability,
#'   the step construction will end \eqn{ \frac{\min(f(x_i), f(x_{i+1}))}{\max(f(x_i), f(x_{i+1}))} < \theta }.
#' Alternatively, if the user wishes to create the steps within certain limits on the
#' x-axis, they can do so by specifying the proposal grid limits using the \code{grid_range} argument.
#'
#' @return
#' A list containing the following elements related to the proposal distribution:
#' \item{grid_data}{A data frame including the created steps information, such as \code{x} (the beginning of each step on the x-axis),
#'  \code{s_upper} (the step height on the y-axis), \code{p_a} (pre-acceptance probability for each step), and \code{s_upper_lower} (a vector used to re-scale the uniform random number when the sample is accepted).}
#' \item{areas}{A vector containing the areas under the left tail bound, the steps in the middle, and the right tail bound.}
#' \item{steps_number}{A scalar representing the number of steps in the proposal.}
#' \item{sampling_probabilities}{A vector containing the areas under the left tail and the combined area of the left tail and middle steps.}
#' \item{unif_scaler}{A scalar representing the inverse probability of sampling from the step part of the proposal, \eqn{\frac{1}{p(lb < x < rb)}}.
#'  Similar to \code{s_upper_lower} in the \code{grid_data} data frame, this value is used to scale the uniform random value when sampling from the steps part of the proposal.}
#' \item{lt_properties}{A vector including 5 values used when sampling under the proposal's left tail using the ARS (Adaptive Rejection Sampling) method.}
#' \item{rt_properties}{A vector including 6 values used when sampling under the proposal's right tail using the ARS method.}
#' \item{alpha}{A scalar representing the uniform step area.}
#' \item{tails_method}{A string representing the tails sampling method, either 'ARS' for Adaptive Rejection Sampling or 'IT' for Inverse Transform.}
#' \item{grid_bounds}{A vector including the left and right bounds of the target density.}
#' \item{dens_func}{The function passed by the user for the target density \code{f}.}
#'
#'
#' @seealso
#' \code{\link{stors}}
#'
#' @examples
#'
#' # Example 1: Building a Grid for Standard Normal Distribution
#' # This example demonstrates constructing a grid for a standard normal distribution
#' # \( f(x) \sim \mathcal{N}(0,1) \),
#' # and shows the optimization table by setting \code{verbose} to \code{TRUE}.
#'
#' # Define the density function, its logarithm,
#' # and its derivative for the standard normal distribution
#' modes_norm = 0
#' f_norm <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x^2) }
#' h_norm <- function(x) { log(f_norm(x)) }
#' h_prime_norm <- function(x) { -x }
#'
#' # Build the proposal grid for the standard normal distribution
#' norm_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_norm,
#'  f = f_norm, h = h_norm, h_prime = h_prime_norm, verbose = TRUE)
#'
#' # Plot the generated grid
#' plot(norm_grid)
#'
#' # Example 2: Grid for a Bimodal Distribution
#' # This example shows how to build a grid for sampling from a bimodal distribution,
#' #combining two normal distributions
#' # \( f(x) = 0.5 \cdot w_1(x) + 0.5 \cdot w_2(x) \),
#' # where \( w_1(x) \sim \mathcal{N}(0, 1) \) and \( w_2(x) \sim \mathcal{N}(4, 1) \).
#'
#' # Define the bimodal density function
#' f_bimodal <- function(x) {
#'   0.5 * (1 / sqrt(2 * pi)) * exp(-(x^2) / 2) + 0.5 * (1 / sqrt(2 * pi)) * exp(-((x - 4)^2) / 2)
#' }
#' modes_bimodal = c(0.00134865, 3.99865)
#'
#' # Build the proposal grid for the bimodal distribution
#' bimodal_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_bimodal, f = f_bimodal)
#'
#' # Print and plot the bimodal grid
#' print(bimodal_grid)
#' plot(bimodal_grid)
#'
#' # Example 3: Grid with 500 Steps for Bimodal Distribution
#' # This example demonstrates constructing a grid with 500 steps,
#' # for the bimodal distribution used in Example 2.
#'
#' bimodal_grid_500 = build_grid(lb = -Inf, rb = Inf, mode = modes_bimodal, f = f_bimodal, steps = 500)
#'
#' # Print and plot the grid with 500 steps
#' print(bimodal_grid_500)
#'
#' @import digest
#' @export
build_grid <- function(lb = -Inf,
                       rb = Inf,
                       modes = NULL,
                       f = NA,
                       h = NULL,
                       h_prime = NULL,
                       steps = NULL,
                       grid_range = NULL,
                       theta = NULL,
                       target_sample_size = 1000,
                       symmetric = NULL,
                       verbose = FALSE) {
  if (!is.function(f)) {
    stop("Error: 'f' density function must be provided.")
  }
  
  if (is.null(h)) {
    h <- function(x) {
      log(f(x))
    }
  }
  
  if (is.null(h_prime)) {
    h_prime <- stors_prime(modes, h)
  }
  
  grid_param <- list(
    target = list(
      density = f,
      log_density = h,
      log_density_prime = h_prime,
      Cumulitive_density = NULL,
      modes = modes,
      modes_count = length(modes),
      between_minima = NULL,
      right_bound = rb,
      left_bound = lb,
      symmetric = symmetric
    ),
    proposal = list(
      grid_range = grid_range,
      tails_method = "ARS",
      steps = steps,
      optimal_step_area = NULL,
      left_steps_proportion = NULL,
      right_steps_proportion = NULL,
      pre_acceptance_threshold = theta,
      target_sample_size = target_sample_size
    ),
    built_in = FALSE,
    cnum = NULL,
    verbose = verbose
  )
  
  
  
  grid_param <- grid_error_checking_and_preparation(grid_param)
  
  grid_param <- grid_check_symmetric(grid_param)
  
  if (!is.null(steps))
    grid_param$proposal$pre_acceptance_thres_hold <- 0.1
  
  optimal_grid_params = find_optimal_grid(grid_param)
  opt_grid <- build_final_grid(gp = optimal_grid_params)
  
  opt_grid$dens_func <- f
  
  class(opt_grid) <- "grid"
  
  if (!(digest(opt_grid) %in% stors_env$created_girds_Id))
    stors_env$created_girds_Id  = append(stors_env$created_girds_Id , digest(opt_grid))
  
  return(opt_grid)
  
}


#' @importFrom utils head
build_final_grid <- function(gp, opt_area = NULL) {
  
  if (is.null(opt_area))
    opt_area <- gp$proposal$optimal_step_area
  
  lb <- gp$target$left_bound
  rb <- gp$target$right_bound
  mode_n <- gp$target$modes_count
  modes <- gp$target$modes
  lstpsp <- gp$proposal$left_steps_proportion
  rstpsp <- gp$proposal$right_steps_proportion
  total_steps <- gp$proposal$steps
  f <- gp$target$density
  h <- gp$target$log_density
  h_prime <- gp$target$log_density_prime
  cdf <- gp$target$Cumulitive_density
  tails_method <- gp$proposal$tails_method
  symmetric <- gp$target$symmetric
  f_params <- gp$f_params
  
  grid_bounds <- rep(NA, 2)
  
  final_grid <- data.frame(
    x = c(),
    s_upper = c(),
    s_lower = c(),
    p_a = c(),
    s_upper_lower = c()
  )
  
  # grids : include lists of steps around each mode
  grids <- list()
  # includes steps number in each grid
  g_len <- c()
  # left_steps[[i]] includes steps data positioned to the left of the i-th mode
  left_steps <- list() 
   # right_steps[[i]] includes steps data positioned to the right of the i-th mode
  right_steps <- list()
  # area under left tail, steps, right tail
  proposal_areas <- c(0, 0, 0)
  
  for (mode_i in  (1:mode_n)) {
    if ((mode_i != 1) || is.null(lstpsp))
    {
      left_steps[[mode_i]] <- find_left_steps(
        gp = gp,
        area = opt_area,
        mode_i = mode_i,
        steps_lim = Inf
      )
      if (!is.null(lstpsp))
        total_steps <- total_steps -  left_steps[[mode_i]]$steps
    }
    if ((mode_i != mode_n) || is.null(rstpsp))
    {
      right_steps[[mode_i]] <- find_right_steps(
        gp = gp,
        area = opt_area,
        mode_i = mode_i,
        steps_lim = Inf
      )
      
      if (!is.null(rstpsp))
        total_steps <- total_steps -  right_steps[[mode_i]]$steps
    }
  }
  
  # after constructing in between modes steps, using the remaining steps
  #(if number of steps is provided to construct tails)
  if (!is.null(rstpsp) || !is.null(lstpsp))
  {
    steps_lim_left = round(lstpsp * total_steps)
    left_steps[[1]] <- find_left_steps(
      gp = gp,
      area = opt_area,
      mode_i = 1,
      steps_lim = steps_lim_left
    )
    steps_lim_right = total_steps - steps_lim_left
    right_steps[[mode_n]] <- find_right_steps(
      gp = gp,
      area = opt_area,
      mode_i = mode_n,
      steps_lim = steps_lim_right
    )
  }
  
  for (mode_i in (1:mode_n)) {
    grids[[mode_i]] <- list(
      data = rbind(left_steps[[mode_i]]$data, right_steps[[mode_i]]$data),
      steps = left_steps[[mode_i]]$steps + right_steps[[mode_i]]$steps
    )
  }
  
  if (mode_n > 1) {
    # binding grids data ( if the minima is not provided)
    for (i in (1:(mode_n - 1))) {
      if (grids[[i]]$data$x[grids[[i]]$steps + 1] > grids[[i + 1]]$data$x[1]) {
        if (grids[[i]]$data$s_upper[grids[[i]]$steps] > grids[[i + 1]]$data$s_upper[1]) {
          grids[[i]]$data <- head(grids[[i]]$data, -1)
          grids[[i]]$data$s_upper[grids[[i]]$steps] <- opt_area / (grids[[i + 1]]$data$x[1] - grids[[i]]$data$x[grids[[i]]$steps])
          grids[[i]]$data$p_a[grids[[i]]$steps] <- 0
          grids[[i + 1]]$data$p_a[1] <- 0
        } else {
          grids[[i + 1]]$data$x[1] <- grids[[i]]$data$x[grids[[i]]$steps + 1]
          grids[[i]]$data <- head(grids[[i]]$data, -1)
          grids[[i + 1]]$data$s_upper[1] <- opt_area / (grids[[i + 1]]$data$x[2] - grids[[i + 1]]$data$x[1])
          grids[[i + 1]]$data$p_a[1] <- 0
          grids[[i]]$data$p_a[grids[[i]]$steps] <- 0
        }
      } else {
        grids[[i]]$steps <- grids[[i]]$steps + 1
        grids[[i]]$data$s_upper[grids[[i]]$steps] <- opt_area / (grids[[i + 1]]$data$x[1] - grids[[i]]$data$x[grids[[i]]$steps])
        grids[[i]]$data$p_a[grids[[i]]$steps] <- 0
      }
      
      m1 <- grids[[i]]$steps
      m2 <- grids[[i + 1]]$steps
      proposal_areas[2] <- proposal_areas[2] + m1 * opt_area
      g_len[length(g_len) + 1] <- m1
      
      if (i == (mode_n - 1)) {
        final_grid <- rbind(final_grid, grids[[i]]$data, grids[[i + 1]]$data)
        proposal_areas[2] <- proposal_areas[2] + m2 * opt_area
        g_len[length(g_len) + 1] <- m2
      } else {
        final_grid <- rbind(final_grid, grids[[i]]$data)
      }
    }
    
  } else {
    final_grid <- grids[[1]]$data
    proposal_areas[2] <- grids[[1]]$steps * opt_area
    g_len[length(g_len) + 1] <- grids[[1]]$steps
  }
  
  steps_number <- sum(g_len)
  x1 <- final_grid$x[1]
  xm <- final_grid$x[steps_number + 1]
  
   if(!is.null(symmetric)) steps_number = steps_number * 2

  if (identical(tails_method, "ARS")) {
    tails_area = tails_ars(final_grid, f, h, h_prime, modes, lb, rb)
    proposal_areas[1] <- tails_area$lta
    proposal_areas[3] <- tails_area$rta
  } else{
    proposal_areas[1] <- cdf(x1)
    proposal_areas[3] <- 1 - cdf(xm)
  }
  
  normalizing_con <- sum(proposal_areas)
  area_cum_sum <- cumsum(proposal_areas)
  sampling_probabilities <- (area_cum_sum / normalizing_con)[1:2]
  unif_scaler <- normalizing_con / proposal_areas[2]
  
  lt_properties <- rep(0, 5)
  rt_properties <- rep(0, 6)
  
  if (identical(tails_method, "ARS")) {
    if (proposal_areas[1] != 0) {
      lt_properties <- c(exp(h_upper(x1, lb, h_prime, h)),
                         normalizing_con * h_prime(x1),
                         h(x1),
                         1 / h_prime(x1),
                         h_prime(x1)) 
    }
    
    if (proposal_areas[3] != 0) {
      rt_properties <- c(
        normalizing_con,
        area_cum_sum[2],
        h_prime(xm) / f(xm),
        1 / h_prime(xm),
        h_prime(xm),
        h(xm)
      )
    }
  }
  
  grid_bounds[1] = lb
  grid_bounds[2] = rb
  
  if(is.null(symmetric)) is_symmetric = FALSE else is_symmetric = TRUE
  
  invisible(
    list(
      grid_data = final_grid,
      areas = proposal_areas,
      steps_number = steps_number,
      sampling_probabilities = sampling_probabilities,
      unif_scaler = unif_scaler,
      lt_properties = lt_properties,
      rt_properties = rt_properties,
      alpha = opt_area,
      tails_method = tails_method,
      grid_bounds = grid_bounds,
      symmetric = symmetric,
      is_symmetric = is_symmetric,
      f_params = f_params
    )
  )
  
}


#' @noRd
find_left_steps <- function(gp, area, mode_i, steps_lim = Inf) {
  lb <- gp$target$left_bound
  rb <- gp$target$right_bound
  mode <- gp$target$modes[mode_i]
  mode_n <- gp$target$modes_count
  f <- gp$target$density
  theta <- gp$proposal$pre_acceptance_threshold
  grid_range <- gp$proposal$grid_range
  # to check if new constructed steps get less than mode_previous due to low density value compared to area
  mode_previous <- ifelse(mode_i == 1, NA, gp$target$modes[mode_i - 1])
  
  memory_res <- (max(500, ceiling(1 / area)) + 500)
  
  x <- rep(NA, memory_res)
  s_upper <- rep(NA, memory_res)
  s_lower <- rep(NA, memory_res)
  p_a <- rep(NA, memory_res)
  s_upper_lower <- rep(NA, memory_res)
  
  l <- 0
  
  i <- memory_res
  
  if (mode != lb) {
    x_previous <- mode
    f_x_previous <- f(x_previous)
    while (TRUE) {
      x_c <- x_previous - area / f_x_previous
      if (l >= steps_lim ||
          x_c < lb ||
          (mode_i == 1 &&
           ((f(x_c) / f_x_previous  <= theta) ||
            x_previous < grid_range[1]))) {
        break
      }
      f_x <- f(x_c)
      if (f(x_c) > f_x_previous ||
          (!is.na(mode_previous) && x_c <  mode_previous)) {
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
  
  d <- data.frame(
    x = x,
    s_upper = s_upper,
    p_a = p_a,
    s_upper_lower = s_upper_lower
  )
  
  d <- subset(d, rowSums(is.na(d)) != ncol(d))
  
  return(list(data = d, steps = l))
}


#' @noRd
find_right_steps <- function(gp,
                             area,
                             mode_i,
                             steps_lim = Inf
                             ) {
  
  lb <- gp$target$left_bound
  rb <- gp$target$right_bound
  mode <- gp$target$modes[mode_i]
  mode_n <- gp$target$modes_count
  f <- gp$target$density
  theta <- gp$proposal$pre_acceptance_threshold
  grid_range <- gp$proposal$grid_range
  
  # to check if x_next exceeds mode_next due to low density value compared to area
  mode_next <- ifelse(mode_i == mode_n, NA, gp$target$modes[mode_i + 1])
  
  memory_res <- (max(500, ceiling(1 / area)) + 500)
  
  x <- rep(NA, memory_res)
  s_upper <- rep(NA, memory_res)
  s_lower <- rep(NA, memory_res)
  p_a <- rep(NA, memory_res)
  s_upper_lower <- rep(NA, memory_res)
  
  r <- 1
  
  if (mode != rb) {
    x_c <- mode
    f_x <- f(x_c)
    while (TRUE) {
      x_next <- x_c + area / f_x
      if (r > steps_lim ||
          x_next > rb ||
          (mode_i == mode_n &&
           ((f(x_next) / f_x <= theta) ||
            x_c > grid_range[2]))) {
        x[r] <- x_c
        s_upper[r] <- s_lower[r] <- s_upper_lower[r] <- p_a[r] <- NA
        break
      }
      f_x_next <- f(x_next)
      if ((f_x_next > f_x) ||
          (!is.na(mode_next) && x_next >  mode_next)) {
        x[r] <- x_c
        s_upper[r] <- f_x
        r_tail_area <- 0
        s_lower[r] <- s_upper_lower[r] <- p_a[r] <- NA
        break
      }
      x[r] <- x_c
      s_upper[r] <- f_x
      s_lower[r] <- f_x_next
      s_upper_lower[r] <- s_upper[r] / s_lower[r]
      p_a[r] <- s_lower[r] / s_upper[r]
      f_x <- f_x_next
      x_c <- x_next
      r <- r + 1
    }
  }
  
  d <- data.frame(
    x = x,
    s_upper = s_upper,
    p_a = p_a,
    s_upper_lower = s_upper_lower
  )
  
  d <- subset(d, rowSums(is.na(d)) != ncol(d))
  
  return(list(data = d, steps = r - 1))
}


#' @noRd
h_upper <- function(grid_point, val, h_prime, h) {
  h_prime(grid_point) * (val - grid_point) + h(grid_point)
}


#' @noRd
tails_ars <- function(grid, f, h, h_prime, modes, lb, rb) {
  steps = length(grid$x)
  modes_n = length(modes)
  # left tail
  if (lb == Inf) {
    l_tail_area <- (1 / h_prime(grid$x[1])) * f(grid$x[1])
  } else{
    if(modes[1] == grid$x[1]){
      l_tail_area <- 0
    }else{
    l_tail_area <- (1 / h_prime(grid$x[1])) * (f(grid$x[1]) - exp(h_upper(grid$x[1], lb, h_prime, h)))
    }
  }
  #right tails
  if (rb == Inf) {
    r_tail_area <- (1 / h_prime(grid$x[steps])) * -f(grid$x[steps])
  } else {
    if(modes[modes_n] == grid$x[steps]){
      r_tail_area <- 0
    }else{
    r_tail_area <- (1 / h_prime(grid$x[steps])) * (exp(h_upper(grid$x[steps], rb, h_prime, h)) - f(grid$x[steps]))
    }
  }
  
  return(list(lta = l_tail_area, rta = r_tail_area))
}


#' @noRd
stors_prime <- function(mode, h) {
  function(x) {
    if (x < mode[1]) {
      x0 = x + 0.0000001
    } else if (x > utils::tail(mode, 1)) {
      x0 = x - 0.0000001
    } else{
      if(x == mode[1])
        warning("First step in the grid starts at mode[1].\nNo Adaptive Rejection Sampling (ARS) for left tail.")
      if(x == utils::tail(mode, 1))
        warning("Last step in the grid ends at mode[n].\nNo ARS for right tail.")
      
      return(0)
    }
    (h(x0) - h(x)) / (x0 - x)
  }
  
}


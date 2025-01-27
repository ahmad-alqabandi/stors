#' Build Proposal
#'
#' Constructs the step optimized proposal density, squeezing function, and log-linear tail proposal for a user defined probability density function.
#'
#' This function is the starting point when a user wishes to build a custom sampler using StORS.
#' It is responsible for generating the step optimized proposal density, squeezing function, and log-linear tail proposal that can be utilized for this purpose.
#' The minimum information that must be supplied by the user is:
#'
#' - The (closed) interval of support for the distribution, \[`lower`, `upper`\] \eqn{\in \mathbb{R}}, which may also be half-closed on either side, or all of \eqn{\mathbb{R}}.
#' - The probability density function (pdf), which need not be normalised, `f`.
#' - Any modes of the pdf, as vector `modes`.
#'
#' Optionally, the log-pdf and derivative of the log-pdf may be supplied.
#'
#' **Arguments for pdf**
#'
#' The pdf (and log-pdf and first derivative of the log-pdf) may depend on certain parameters.
#' If so, these can be from the second argument onward in `f`.
#' For instance, consider the Kumaraswamy distribution, which has pdf:
#' \deqn{f(x; a,b) = a b x^{a-1}{ (1-x^a)}^{b-1},  \ \ \mbox{where} \ \ x \in (0,1)}
#' This pdf has known modes.
#'
#' Then, to implement as a custom StORS sampler, we would first define the pdf in R:
#'
#' \code{dkumaraswamy <- function(x, a, b) a*b*(x^(a-1))*(1-x^a)^(b-1)}
#'
#' Then, to construct a StORS proposal for \eqn{a=2} and \eqn{b=2}, we would call
#'
#' \code{Proposal <- build_Proposal(lower = 0, upper = 1, modes = sqrt(1/3), f = dkumaraswamy, a = 2, b = 2)}
#'
#' **StORS proposal construction**
#'
#' StORS defines an unnormalised piecewise constant proposal density and squeezing function, with a Proposal defining the change points.
#' To optimise the execution speed on modern CPUs, the unnormalised piecewise constant proposal has fixed area for each segment with one end of the segment coinciding with the user's pdf.
#' That is, each step of the function has width defined by \eqn{w_i = (x_i - x_{i-1})} and a height determined by \eqn{h_i = \max(f(x_{i-1}), f(x_i))}, such that \eqn{w_i h_i = \alpha \ \forall\,i} where \eqn{\alpha} is constant.
#'
#' Once the user has constructed the proposal, the sampling function can be built using [build_sampler()].
#'
#' **Internal details**
#'
#' The function \code{build_final_Proposal()} manages the construction of these steps and calculates values critical for the sampling process.
#' When the resultant Proposal is used with the \code{build_sampler()} function, these values are cached,
#' significantly enhancing the computational efficiency and hence improving sampling speed.
#'  During the optimization process, we aim for a certain Proposal
#' size based on L1-3 memory cache size. Therefore, we test the speed of Proposals of sizes \eqn{2^m} Kb.
#'  To achieve this, we estimate the uniform step area
#' based on a certain steps number that leads to the target cache size,
#'  \eqn{ \alpha = \frac{1}{\text{number of steps}} }.
#'
#'
#' The speed testing for each possible Proposal is initially based on a sample size of 1000.
#'  However, if the user wishes to optimize the Proposal for a different sample size, they can do so
#' by specifying the desired sample size using the \code{target_sample_size} argument.
#'
#' In case the user wants to select a specific number of steps for the proposal
#' and bypass the optimization process, this can be done by specifying a steps number greater than the number of modes by 2 using the \code{steps} argument.
#'  If the target density is heavy-tailed,
#'   and the user wishes to stop the Proposal building process at a certain pre-acceptance threshold, this can be achieved by setting
#' the acceptance probability threshold \code{theta} \eqn{\theta}.
#'  Once the steps reach this level of pre-acceptance probability,
#'   the step construction will end \eqn{ \frac{\min(f(x_i), f(x_{i+1}))}{\max(f(x_i), f(x_{i+1}))} < \theta }.
#' Alternatively, if the user wishes to create the steps within certain limits on the
#' x-axis, they can do so by specifying the proposal limits using the \code{proposal_range} argument.
#'
#'#' @param modes
#'        Numeric vector of modes of the density function.
#' @param f
#'        A function which returns the (unnormalised) probability density function of the target distribution.
#'        The first argument must be the value at which the pdf is to be evaluated.
#'        Additional arguments may be parameters of the distribution, which should be specified by name in the `...` arguments.
#' @param lower
#'        Numeric scalar representing the lower bound of the target density.
#'        Default is `-Inf` for unbounded lower support.
#' @param upper
#'        Numeric scalar representing the upper bound of the target density.
#'        Default is `Inf` for unbounded upper support.
#' @param h
#'        An optional function which returns the (unnormalised) log-probability density function of the target distribution.
#'        As for `f` the first argument must be the value at which the log-pdf is to be evaluated and additional parameters may be named arguments passed to `...`.
#' @param h_prime
#'        An optional function which returns the first derivative of the (unnormalised) log-probability density function of the target distribution.
#'        As for `f` the first argument must be the value at which the log-pdf is to be evaluated and additional parameters may be named arguments passed to `...`.
#' @param steps
#'        Optional integer scalar specifying the number of steps in the step optimised part of the proposal density and squeezing function.
#' @param proposal_range
#'        Optional numeric vector of length 2 specifying the lower and upper range of the steps in the step optimised part of the proposal density and squeezing function.
#'        This range should be contained within the interval defined by `lower` and `upper`.
#' @param theta
#'        Optional numeric scalar (between 0 and 1) defining the pre-acceptance threshold.
#'        This dictates when no further steps should be added in the step optimised part of the proposal density and squeezing function, based on the probability of pre-acceptance.
#' @param target_sample_size
#'        Integer scalar indicating the typical sample size that will be requested when sampling from this density using build_sampler.
#'        The proposal optimization process bases benchmark timings on this target size in order to select a proposal best suited to the desired sample size.
#'        Note this does *not* limit sampling to this number, it is merely a guide should the user be aware that a certain sample size will be most commonly sampled.
#' @param verbose
#'        Logical scalar.
#'        If `TRUE`, a table detailing the optimization areas and steps will be displayed during proposal optimization.
#'        Defaults to `FALSE`.
#' @param ...
#'        Further arguments to be passed to `f`, `h`, and `h_prime`, if they depend on additional parameters.
#'
#' @return
#' This returns a list which is used to construct the sampler by passing to \code{\link{build_sampler}} function.
#'
#' A list containing the optimized proposal and related parameters for the specified built-in distribution:
#' \describe{
#'   \item{\code{data}}{A data frame with detailed information about the proposal steps, including:
#'   \describe{
#'     \item{\code{x}}{The start point of each step on the x-axis.}
#'     \item{\code{s_upper}}{The height of each step on the y-axis.}
#'     \item{\code{p_a}}{Pre-acceptance probability for each step.}
#'     \item{\code{s_upper_lower}}{A vector used to scale the uniform random number when the sample is accepted.}
#'   }}
#'   \item{\code{areas}}{A numeric vector containing the areas under:
#'   \describe{
#'     \item{\code{left_tail}}{The left tail bound.}
#'     \item{\code{steps}}{The middle steps.}
#'     \item{\code{right_tail}}{The right tail bound.}
#'   }}
#'   \item{\code{steps_number}}{An integer specifying the number of steps in the proposal.}
#'   \item{\code{sampling_probabilities}}{A numeric vector with:
#'   \describe{
#'     \item{\code{left_tail}}{The probability of sampling from the left tail.}
#'     \item{\code{left_and_middle}}{The combined probability of sampling from the left tail and middle steps.}
#'   }}
#'   \item{\code{unif_scaler}}{A numeric scalar, the inverse probability of sampling from the steps part of the proposal (\eqn{\frac{1}{p(lower < x < upper)}}). Used for scaling uniform random values.}
#'   \item{\code{lt_properties}}{A numeric vector of 5 values required for Adaptive Rejection Sampling (ARS) in the left tail.}
#'   \item{\code{rt_properties}}{A numeric vector of 6 values required for ARS in the right tail.}
#'   \item{\code{alpha}}{A numeric scalar representing the uniform step area.}
#'   \item{\code{tails_method}}{A string, either \code{"ARS"} (Adaptive Rejection Sampling) or \code{"IT"} (Inverse Transform), indicating the sampling method for the tails.}
#'   \item{\code{proposal_bounds}}{A numeric vector specifying the left and right bounds of the target density.}
#'   \item{\code{cnum}}{An integer representing the cache number of the created proposal in memory.}
#'   \item{\code{symmetric}}{A numeric scalar indicating the symmetry point of the proposal, or \code{NULL} if not symmetric.}
#'   \item{\code{f_params}}{A list of parameters for the target density that the proposal is designed for.}
#'   \item{\code{is_symmetric}}{A logical value indicating whether the proposal is symmetric.}
#'   \item{\code{proposal_type}}{A string indicating the type of the generated proposal:
#'   \describe{
#'     \item{\code{"scaled"}}{The proposal is "scalable" and standardized with \code{rate = 1}. This is used when parameter \code{rate} is either \code{NULL} or not provided. Scalable proposals are compatible with \code{\link{srexp}}.}
#'     \item{\code{"custom"}}{The proposal is "custom" when \code{rate} is provided. Custom proposals are compatible with \code{\link{srexp_custom}}.}
#'   }}
#'   \item{\code{target_function_area}}{A numeric scalar estimating the area of the target distribution.}
#'   \item{\code{dens_func}}{A string containing the hardcoded density function.}
#'   \item{\code{density_name}}{A string specifying the name of the target density distribution.}
#'   \item{\code{lock}}{An identifier used for saving and loading the proposal from disk.}
#' }
#'
#'
#' @seealso
#' \code{\link{build_sampler}}: Function to build and return a sampling function based on the provided proposal properties.
#'
#'
#' @examples
#'
#' # Example 1: Building a proposal for Standard Normal Distribution
#' # This example demonstrates constructing a proposal for a standard normal distribution
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
#' # Build the proposal for the standard normal distribution
#' norm_proposal = build_proposal(lower = -Inf, upper = Inf, mode = modes_norm,
#'  f = f_norm, h = h_norm, h_prime = h_prime_norm, verbose = TRUE)
#'
#' # Plot the generated proposal
#' plot(norm_proposal)
#'
#' # Example 2: proposal for a Bimodal Distribution
#' # This example shows how to build a proposal for sampling from a bimodal distribution,
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
#' # Build the proposal for the bimodal distribution
#' bimodal_proposal = build_proposal(lower = -Inf, upper = Inf, mode = modes_bimodal, f = f_bimodal)
#'
#' # Print and plot the bimodal proposal
#' print(bimodal_proposal)
#' plot(bimodal_proposal)
#'
#' # Example 3: Proposal with 500 Steps for Bimodal Distribution
#' # This example demonstrates constructing a proposal with 500 steps,
#' # for the bimodal distribution used in Example 2.
#'
#' bimodal_proposal_500 = build_proposal(lower = -Inf, upper = Inf, mode = modes_bimodal, f = f_bimodal, steps = 500)
#'
#' # Print and plot the proposal with 500 steps
#' print(bimodal_proposal_500)
#'
#' @import digest
#' @export
build_proposal <- function(f = NULL,
                           modes = NA,
                           lower = -Inf,
                           upper = Inf,
                           h = NULL,
                           h_prime = NULL,
                           steps = NULL,
                           proposal_range = NULL,
                           theta = NULL,
                           target_sample_size = 1000,
                           verbose = FALSE,
                           ...) {
  density_arguments <- list(...)


  if (!is.function(f)) {
    stop("Error: 'f' density function must be provided.")
  }

  preserved_f <- f
  f <- create_function(f, density_arguments)

  if (is.null(h)) {
    h <- function(x) {
      log(f(x))
    }
  } else {
    h <- create_function(h, density_arguments)
  }

  if (is.null(h_prime)) {
    h_prime <- estimate_slope(modes, h)
  } else {
    h_prime <- create_function(h_prime, density_arguments)
  }


  modes <- adjust_modes(modes, lower, upper, f)

  proposal_param <- list(
    target = list(
      density = f,
      density_arguments = density_arguments,
      log_density = h,
      log_density_prime = h_prime,
      Cumulitive_density = NULL,
      modes = modes,
      modes_count = length(modes),
      between_minima = NULL,
      right_bound = upper,
      left_bound = lower,
      estimated_area = NULL,
      symmetric = NULL
    ),
    proposal = list(
      proposal_range = proposal_range,
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
    proposal_type = NULL,
    c_function_name = NULL,
    verbose = verbose,
    f_params = NULL
  )



  proposal_param <- proposal_error_checking_and_preparation(proposal_param)

  if (!is.null(steps))
    proposal_param$proposal$pre_acceptance_thres_hold <- 0.1

  optimal_proposal_params <- find_optimal_proposal(proposal_param)
  opt_proposal <- build_final_proposal(gp = optimal_proposal_params)

  opt_proposal$dens_func <- deparse(preserved_f)

  lock <- digest(opt_proposal)
  opt_proposal$lock <- lock

  class(opt_proposal) <- "proposal"

  return(opt_proposal)

}




construct_left_and_right_steps <- function(gp, opt_area, mode_n) {
  total_steps <- gp$proposal$steps
  lstpsp <- gp$proposal$left_steps_proportion
  rstpsp <- gp$proposal$right_steps_proportion
  left_steps <- list()
  right_steps <- list()
  grids <- list()

  for (mode_i in  (1:mode_n)) {
    if ((mode_i != 1) || is.null(lstpsp)) {
      left_steps[[mode_i]] <- find_left_steps(
        gp = gp,
        area = opt_area,
        mode_i = mode_i,
        steps_lim = Inf
      )
      if (!is.null(lstpsp))
        total_steps <- total_steps -  left_steps[[mode_i]]$steps
    }
    if ((mode_i != mode_n) || is.null(rstpsp)) {
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

  if (!is.null(rstpsp) || !is.null(lstpsp)) {
    steps_lim_left <- round(lstpsp * total_steps)
    left_steps[[1]] <- find_left_steps(
      gp = gp,
      area = opt_area,
      mode_i = 1,
      steps_lim = steps_lim_left
    )
    steps_lim_right <- total_steps - steps_lim_left
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
  return(grids)

}


#' @importFrom utils head
build_final_proposal <- function(gp, opt_area = NULL) {
  if (is.null(opt_area))
    opt_area <- gp$proposal$optimal_step_area

  lower <- gp$target$left_bound
  upper <- gp$target$right_bound
  mode_n <- gp$target$modes_count
  modes <- gp$target$modes
  f <- gp$target$density
  h <- gp$target$log_density
  h_prime <- gp$target$log_density_prime
  cdf <- gp$target$Cumulitive_density
  tails_method <- gp$proposal$tails_method
  symmetric <- gp$target$symmetric
  f_params <- gp$f_params
  proposal_type <- gp$proposal_type
  cnum <- gp$cnum
  f_area <- gp$target$estimated_area
  density_arguments <- gp$target$density_arguments



  proposal_bounds <- rep(NA, 2)

  final_proposal <- data.frame(
    x = c(),
    s_upper = c(),
    s_lower = c(),
    p_a = c(),
    s_upper_lower = c()
  )

  g_len <- c()
  proposal_areas <- c(0, 0, 0)

  grids <- construct_left_and_right_steps(gp, opt_area, mode_n)

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
        final_proposal <- rbind(final_proposal, grids[[i]]$data, grids[[i + 1]]$data)
        proposal_areas[2] <- proposal_areas[2] + m2 * opt_area
        g_len[length(g_len) + 1] <- m2
      } else {
        final_proposal <- rbind(final_proposal, grids[[i]]$data)
      }
    }

  } else {
    final_proposal <- grids[[1]]$data
    proposal_areas[2] <- grids[[1]]$steps * opt_area
    g_len[length(g_len) + 1] <- grids[[1]]$steps
  }

  steps_number <- sum(g_len)
  x1 <- final_proposal$x[1]
  xm <- final_proposal$x[steps_number + 1]


  if (!is.null(symmetric))
    steps_number <- steps_number * 2

  if (identical(tails_method, "ARS")) {
    tails_area <- tails_ars(final_proposal, f, h, h_prime, modes, lower, upper)

    if (modes[1] == lower)
      proposal_areas[1] <- 0
    else
      proposal_areas[1] <- tails_area$lta

    if (modes[mode_n] == upper)
      proposal_areas[3] <- 0
    else
      proposal_areas[3] <- tails_area$rta

  } else {
    if (modes[1] == lower) {
      proposal_areas[1] <- 0
    } else {
      proposal_areas[1] <- cdf(x1) - cdf(lower)
    }

    if (modes[mode_n] == upper) {
      proposal_areas[3] <- 0
    } else {
      proposal_areas[3] <- cdf(upper) - cdf(xm)
    }
  }

  normalizing_con <- sum(proposal_areas)
  area_cum_sum <- cumsum(proposal_areas)
  sampling_probabilities <- (area_cum_sum / normalizing_con)[1:2]
  unif_scaler <- normalizing_con / proposal_areas[2]

  lt_properties <- rep(0, 5)
  rt_properties <- rep(0, 6)

  if (identical(tails_method, "ARS")) {
    if (proposal_areas[1] != 0) {
      lt_properties <- c(
        exp(h_upper(x1, lower, h_prime, h)),
        normalizing_con * h_prime(x1),
        h(x1),
        1 / h_prime(x1),
        h_prime(x1)
      )
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

  proposal_bounds[1] <- lower
  proposal_bounds[2] <- upper

  if (is.null(symmetric))
    is_symmetric <- FALSE
  else
    is_symmetric <- TRUE

  # if(steps_number != length(final_proposal$x) - 1 ) stop("ERROR: STEPS NUMBER")

  invisible(
    list(
      data = final_proposal,
      areas = proposal_areas,
      steps_number = steps_number,
      sampling_probabilities = sampling_probabilities,
      unif_scaler = unif_scaler,
      lt_properties = lt_properties,
      rt_properties = rt_properties,
      alpha = opt_area,
      tails_method = tails_method,
      proposal_bounds = proposal_bounds,
      cnum = cnum,
      symmetric = symmetric,
      f_params = f_params,
      is_symmetric = is_symmetric,
      proposal_type = proposal_type,
      target_function_area = f_area,
      density_arguments = density_arguments
    )
  )

}


#' @noRd
find_left_steps <- function(gp, area, mode_i, steps_lim = Inf) {
  estimated_area <- gp$target$estimated_area
  lower <- gp$target$left_bound
  mode <- gp$target$modes[mode_i]
  f <- gp$target$density
  theta <- gp$proposal$pre_acceptance_threshold
  proposal_range <- gp$proposal$proposal_range
  # to check if new constructed steps get less than mode_previous due to low density value compared to area
  mode_previous <- ifelse(mode_i == 1, NA, gp$target$modes[mode_i - 1])

  init_memory_res <- (min(1000, ceiling(estimated_area / area)) + 500) * 2

  x <- rep(NA, init_memory_res)
  s_upper <- rep(NA, init_memory_res)
  s_lower <- rep(NA, init_memory_res)
  p_a <- rep(NA, init_memory_res)
  s_upper_lower <- rep(NA, init_memory_res)

  i <- 0

  if (mode != lower) {
    x_previous <- mode
    f_x_previous <- f(x_previous)
    while (TRUE) {
      x_c <- x_previous - area / f_x_previous

      if (i >= steps_lim ||
          x_c < lower ||
          (mode_i == 1 &&
           ((f(x_c) / f_x_previous  <= theta) ||
            x_previous < proposal_range[1])))
        break

      f_x <- f(x_c)

      if (f(x_c) > f_x_previous ||
          (!is.na(mode_previous) && x_c <  mode_previous))
        break

      i <- i + 1
      x[i] <- x_c
      s_upper[i] <- f_x_previous
      s_lower[i] <- f_x
      s_upper_lower[i] <- s_upper[i] / s_lower[i]
      p_a[i] <- f_x / f_x_previous
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

  d <- d[rev(seq_along(row(d))), ]

  d <- subset(d, rowSums(is.na(d)) != ncol(d))

  return(list(data = d, steps = i))
}


#' @noRd
find_right_steps <- function(gp, area, mode_i, steps_lim = Inf) {
  estimated_area <- gp$target$estimated_area
  upper <- gp$target$right_bound
  mode <- gp$target$modes[mode_i]
  mode_n <- gp$target$modes_count
  f <- gp$target$density
  theta <- gp$proposal$pre_acceptance_threshold
  proposal_range <- gp$proposal$proposal_range

  # to check if x_next exceeds mode_next due to low density value compared to area
  mode_next <- ifelse(mode_i == mode_n, NA, gp$target$modes[mode_i + 1])

  init_memory_res <- (min(1000, ceiling(estimated_area / area)) + 500)

  x <- rep(NA, init_memory_res)
  s_upper <- rep(NA, init_memory_res)
  s_lower <- rep(NA, init_memory_res)
  p_a <- rep(NA, init_memory_res)
  s_upper_lower <- rep(NA, init_memory_res)

  r <- 1

  if (mode != upper) {
    x_c <- mode
    f_x <- f(x_c)
    while (TRUE) {
      x_next <- x_c + area / f_x
      if (r > steps_lim ||
          x_next > upper ||
          (mode_i == mode_n &&
           ((f(x_next) / f_x <= theta) ||
            x_c > proposal_range[2]))) {
        x[r] <- x_c
        s_upper[r] <- s_lower[r] <- s_upper_lower[r] <- p_a[r] <- NA
        break
      }
      f_x_next <- f(x_next)
      if ((f_x_next > f_x) ||
          (!is.na(mode_next) && x_next >  mode_next)) {
        x[r] <- x_c
        s_upper[r] <- f_x
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
  } else {
    x[r] <- mode
    s_upper[r] <- s_lower[r] <- s_upper_lower[r] <-  p_a[r] <- 0
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
h_upper <- function(proposals_point, val, h_prime, h) {
  h_prime(proposals_point) * (val - proposals_point) + h(proposals_point)
}


#' @noRd
tails_ars <- function(grid, f, h, h_prime, modes, lower, upper) {
  steps <- length(grid$x)
  modes_n <- length(modes)
  # left tail
  if (lower == Inf) {
    l_tail_area <- (1 / h_prime(grid$x[1])) * f(grid$x[1])
  } else {
    if (modes[1] == grid$x[1]) {
      l_tail_area <- 0
    } else {
      l_tail_area <- (1 / h_prime(grid$x[1])) * (f(grid$x[1]) - exp(h_upper(grid$x[1], lower, h_prime, h)))
    }
  }
  #right tails
  if (upper == Inf) {
    r_tail_area <- (1 / h_prime(grid$x[steps])) * -f(grid$x[steps])
  } else {
    if (modes[modes_n] == grid$x[steps]) {
      r_tail_area <- 0
    } else {
      r_tail_area <- (1 / h_prime(grid$x[steps])) * (exp(h_upper(grid$x[steps], upper, h_prime, h)) - f(grid$x[steps]))
    }
  }

  return(list(lta = l_tail_area, rta = r_tail_area))
}


#' @noRd
estimate_slope <- function(mode, h) {
  function(x) {
    if (x < mode[1]) {
      x0 <- x + 0.0000001
    } else if (x > utils::tail(mode, 1)) {
      x0 <- x - 0.0000001
    } else {
      if (x == mode[1])
        warning(
          "First step in the grid starts at mode[1].\nNo Adaptive Rejection Sampling (ARS) for left tail."
        )
      if (x == utils::tail(mode, 1))
        warning("Last step in the grid ends at mode[n].\nNo ARS for right tail.")

      return(0)
    }
    (h(x0) - h(x)) / (x0 - x)
  }

}

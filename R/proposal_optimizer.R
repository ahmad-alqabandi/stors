#' @noRd
proposal_optimizer <- function(dendata,
                               density_name,
                               xl = NULL,
                               xr = NULL,
                               f,
                               modes,
                               f_params = NULL,
                               steps = NULL,
                               proposal_range = NULL,
                               theta = 0.1,
                               target_sample_size = NULL,
                               proposal_type,
                               symmetric = NULL,
                               cnum = NULL,
                               verbose = FALSE) {
  free_cache_cnum_c(cnum)

  if ((is.null(xl) || (xl < dendata$lower)))
    xl <- dendata$lower

  if (is.null(xr) || (xr > dendata$upper))
    xr <- dendata$upper

  if (xl > xr) {
    cli::cli_abort(c("x" = "{.strong xl must be smaller than xr.}",
                     "i" = "You provided: xl = {xl}, xr = {xr}."))
  }

  modes <- adjust_modes(modes, xl, xr, f)


  proposal_param <- list(
    target = list(
      density = f,
      log_density = NULL,
      log_density_prime = NULL,
      Cumulitive_density = NULL,
      modes = modes,
      modes_count = length(modes),
      between_minima = NULL,
      right_bound = xr,
      left_bound = xl,
      estimated_area = NULL,
      symmetric = symmetric,
      density_arguments = list()
    ),

    proposal = list(
      proposal_range = proposal_range,
      tails_method = dendata$tails_method,
      steps = steps,
      optimal_step_area = NULL,
      left_steps_proportion = NULL,
      right_steps_proportion = NULL,
      pre_acceptance_threshold = theta,
      target_sample_size = target_sample_size
    ),

    built_in = TRUE,
    cnum = cnum,
    proposal_type = proposal_type,
    c_function_name = density_name,
    verbose = verbose,
    f_params = f_params
  )

  if (identical(dendata$tails_method, "ARS")) {
    h <- function(x) {
      log(f(x))
    }

    h_prime <- estimate_slope(modes, h)

  } else {
    cdf <- do.call(dendata$create_cdf, f_params)

  }

  if (proposal_param$proposal$tails_method == "IT") {
    proposal_param$target$Cumulitive_density <- cdf

  } else {
    proposal_param$target$log_density <- h
    proposal_param$target$log_density_prime <- h_prime

  }

  proposal_param <- proposal_error_checking_and_preparation(proposal_param)
  proposal_param <- proposal_check_symmetric(proposal_param)

  optimal_proposal_params <- find_optimal_proposal(proposal_param)
  opt_proposal <- build_final_proposal(gp = optimal_proposal_params)

  cache_proposal_c(cnum, opt_proposal)
  opt_proposal$dens_func <- deparse(f)
  opt_proposal$density_name <- density_name


  lock <- digest(opt_proposal)
  opt_proposal$lock <- lock

  class(opt_proposal) <- "proposal"

  save_builtin_proposal(cnum, opt_proposal)

  return(opt_proposal)

}




#' @importFrom microbenchmark microbenchmark
#' @importFrom stats integrate
find_optimal_proposal <- function(gp) {
  target_sample_size <- gp$proposal$target_sample_size
  theta <- gp$proposal$pre_acceptance_threshold
  proposal_range <- gp$proposal$proposal_range
  lower <- gp$target$left_bound
  upper <- gp$target$right_bound
  f <- gp$target$density
  lstpsp <- gp$proposal$left_steps_proportion
  rstpsp <- gp$proposal$right_steps_proportion
  mode_n <- gp$target$modes_count
  steps <- gp$proposal$steps
  c_function_name <- gp$c_function_name
  verbose <- gp$verbose
  opt_area <- NULL


  opt_alpha_length <- 3
  opt_times <- 10000
  opt_cache_sizes <- c(4, 8, 16, 32, 64, 128, 256, 512, 1024)
  opt_df_var <- 32 # 4 doubles
  opt_list_var <- 160 # 20 doubles
  opt_steps <- round(((opt_cache_sizes * 1024) - opt_list_var) / opt_df_var)
  opt_areas <- 1 / opt_steps

  times <- ceiling(opt_times / target_sample_size)

  f_integrate <- integrate(f, lower = lower, upper = upper)
  relative_error <- f_integrate$abs.error / f_integrate$value * 100

  if (relative_error > 0.1) {
    cli::cli_abort(c("x" = "The provided density function has a large relative integration error: {sprintf('%.3f', relative_error)} ",
                     "i" = "This may indicate numerical instability or an issue with the function's behavior over the truncation range.",
                     "i" = "Integration was performed over the range: [{lower}, {upper}]. Consider adjusting the truncation limits or reviewing the density function."))
  }


  f_area <- f_integrate$value
  gp$target$estimated_area <- f_area

  if ((theta == 0 &&
       identical(proposal_range, c(lower, upper))) ||
      !is.null(steps)) {
    if (!is.null(steps)) {
      opt_area <- 1 / steps * f_area
    } else {
      opt_area <- 1 / 4096 * f_area
    }

    gdl <- find_left_steps(
      gp = gp,
      area = opt_area,
      mode_i = 1,
      steps_lim = Inf
    )

    max_left_stps <- gdl$steps

    gdr <- find_right_steps(
      gp = gp,
      area = opt_area,
      mode_i = mode_n,
      steps_lim = Inf
    )


    max_right_stps <- gdr$steps


    lstpsp <- max_left_stps / (max_left_stps + max_right_stps)
    rstpsp <- 1 - lstpsp

    if (!is.null(steps)) {
      gp$proposal$optimal_step_area <- opt_area
      gp$proposal$left_steps_proportion <- lstpsp
      gp$proposal$right_steps_proportion <- rstpsp
      return(gp)
    }

  }

  performance <- data.frame(area = numeric(),
                            time = numeric(),
                            steps = numeric())


  if (gp$built_in) {
    cnum <- gp$cnum
    density_fun <- get_buildin_sampling_function(cnum, c_function_name)
  } else {
    cnum <- 0
    rfunc_env <- new.env()
    density_fun <- function(n) {
      .Call(C_stors, n, cnum, f, rfunc_env)
    }
  }

  custom_opt_areas <- opt_areas * f_area

  for (i in seq_along(custom_opt_areas)) {
    area <- custom_opt_areas[i]

    if (is.null(rstpsp)) {
      step <- Inf
    } else {
      step <- opt_steps[i]
    }

    if (verbose) {

      if (i == 1)
        cli::cli_h1("Optimization Summary")

      cat("\n\n",
        sprintf("Step: %10s | Area: %10.9f", step, area),
        "\n-----------------------------------------------\n",
        "       Steps      |       Area       |    Best Sim Time\n",
        "-----------------------------------------------\n"
      )
    }

    area_seq <- seq(from = area * 0.95,
                    to = area * 1.05,
                    length.out = opt_alpha_length)

    steps_time <- double()
    steps_number <- double()

    for (j in seq_along(area_seq)) {
      grid <- build_final_proposal(gp = gp, opt_area = area_seq[j])
      cache_proposal_c(cnum, grid)
      gc <- gc()

      if (gp$built_in) {
        suppressWarnings({
          cost <- microbenchmark::microbenchmark(st = density_fun(target_sample_size),
                                                 times = times)
        })
      } else {
        suppressWarnings({
          cost <- microbenchmark::microbenchmark(st = density_fun(target_sample_size),
                                                 times = times)
        })
      }

      free_cache_cnum_c(cnum)

      steps_time <- append(steps_time, stats::median(cost$time[cost$expr == "st"]))
      steps_number <- append(steps_number, grid$steps_number)

      if (verbose) {
        cat(sprintf("--- %10s ---  --- %10.9f ---   --- %10.2f ---\n", grid$steps_number, area_seq[j], steps_time[j]))
      }
    }

    if (i != 1 &&
        min(steps_time) >= min(performance$time, na.rm = TRUE)) {

      performance <- stats::na.omit(performance)

      if (verbose) {
        cat("\n")
        cli::cli_rule()
        cli::cli_h1("Performance Data:")
        cat("     Area       |     Time     |   Steps\n")
        cat("-----------------------------------------\n")
        for (k in seq_len(nrow(performance))) {
          cat(
            sprintf(
              "%13.10f | %10.2f | %7s\n",
              performance$area[k],
              performance$time[k],
              performance$steps[k]
            )
          )
        }
        cli::cli_rule()
      }
      break
    }

    min_ind <- which(steps_time == min(steps_time, na.rm = TRUE))[1]

    performance[nrow(performance) + 1, ] <- c(area_seq[min_ind], steps_time[min_ind], steps_number[min_ind])

  }

  opt_performance <- performance[which(performance$time == min(performance$time, na.rm = TRUE))[1], ]
  gp$proposal$optimal_step_area <- opt_performance$area
  gp$proposal$steps <- opt_performance$steps
  gp$proposal$left_steps_proportion <- lstpsp
  gp$proposal$right_steps_proportion <- rstpsp

  return(gp)

}

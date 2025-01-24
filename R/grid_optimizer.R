#' @noRd
grid_optimizer <- function(dendata,
                           density_name,
                           xl = NULL,
                           xr = NULL,
                           f,
                           modes,
                           f_params = NULL,
                           steps = NULL,
                           grid_range = NULL,
                           theta = NULL,
                           target_sample_size = NULL,
                           grid_type,
                           symmetric = NULL,
                           cnum = NULL,
                           verbose = FALSE) {

  free_cache_cnum_c(cnum)

  if ((is.null(xl) || (xl < dendata$lb)))
    xl <- dendata$lb

  if (is.null(xr) || (xr > dendata$rb))
    xr <- dendata$rb

  if (xl > xr)
    stop("xl must be smaller than xr.")

  modes <- adjust_modes(modes, xl, xr, f)


  grid_param <- list(
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
      grid_range = grid_range,
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
    grid_type = grid_type,
    c_function_name = density_name,
    verbose = verbose,
    f_params = f_params
  )

  if (identical(dendata$tails_method, "ARS")) {

    h <- function(x) {
      log(f(x))
      }

    h_prime <- stors_prime(modes, h)

  } else {
    cdf <- do.call(dendata$create_cdf, f_params)

  }

  if (grid_param$proposal$tails_method == "IT") {
    grid_param$target$Cumulitive_density <- cdf

  } else {
    grid_param$target$log_density <- h
    grid_param$target$log_density_prime <- h_prime

  }

  grid_param <- grid_error_checking_and_preparation(grid_param)
  grid_param <- grid_check_symmetric(grid_param)

  optimal_grid_params <- find_optimal_grid(grid_param)
  opt_grid <- build_final_grid(gp = optimal_grid_params)

  cache_grid_c(cnum, opt_grid)
  opt_grid$dens_func <- deparse(f)
  opt_grid$density_name <- density_name


  lock <- digest(opt_grid)
  opt_grid$lock <- lock

  class(opt_grid) <- "grid"

  save_builtin_grid(cnum, opt_grid)

  return(opt_grid)

}




#' @importFrom microbenchmark microbenchmark
#' @importFrom stats integrate
find_optimal_grid <- function(gp) {
  target_sample_size <- gp$proposal$target_sample_size
  theta <- gp$proposal$pre_acceptance_threshold
  grid_range <- gp$proposal$grid_range
  lb <- gp$target$left_bound
  rb <- gp$target$right_bound
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

  f_integrate <- integrate(f, lower = lb, upper = rb)
  relative_error <- f_integrate$abs.error / f_integrate$value * 100

  if (relative_error > 0.1) {
    stop(paste0(
      "provided density has large relative error = ",
      relative_error,
      "%"
    ))
  }

  f_area <- f_integrate$value
  gp$target$estimated_area <- f_area

  if ((theta == 0 &&
       identical(grid_range, c(lb, rb))) ||
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

  performance <- data.frame(area = numeric(), time = numeric(), steps = numeric())


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
      cat(
        "\n=====================================\n",
        sprintf("Step: %10s | Area: %10.9f", step, area),
        "\n-------------------------------------\n",
        "       Area       |    Best Sim Time\n",
        "-------------------------------------\n"
      )
    }

    area_seq <- seq(from = area * 0.95,
                   to = area * 1.05,
                   length.out = opt_alpha_length)

    steps_time <- double()

    for (j in seq_along(area_seq)) {
      grid <- build_final_grid(gp = gp, opt_area = area_seq[j])
      cache_grid_c(cnum, grid)
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

      if (verbose) {
        cat(sprintf("--- %10.9f ---   --- %10.2f ---\n", area_seq[j], steps_time[j]))
      }
    }

    if (i != 1 &&
        min(steps_time) >= min(performance$time, na.rm = TRUE)) {
      if (verbose) {
        cat("\n=====================================\n")
        cat("\nPerformance Data:\n")
        cat("     Area       |     Time     |   Steps\n")
        cat("-----------------------------------------\n")
        for (k in seq_along(row(performance))) {
          cat(
            sprintf(
              "%13.10f | %10.2f | %7.2f\n",
              performance$area[k],
              performance$time[k],
              performance$steps[k]
            )
          )
        }
        cat("=====================================\n\n")
      }
      break
    }

    min_ind <- which(steps_time == min(steps_time, na.rm = TRUE))[1]
    performance[nrow(performance) + 1, ] <- c(area_seq[min_ind], steps_time[min_ind], step)

  }

  opt_performance <- performance[which(performance$time == min(performance$time, na.rm = TRUE))[1], ]
  gp$proposal$optimal_step_area <- opt_performance$area
  gp$proposal$steps <- opt_performance$steps
  gp$proposal$left_steps_proportion <- lstpsp
  gp$proposal$right_steps_proportion <- rstpsp

  return(gp)

}

#' @noRd
adjust_modes <- function(mode, xl, xr, f) {
  delta <- 0.0001
  new_modes <- mode[mode >= xl & mode <= xr]

  if (length(new_modes) == 0) {
    if (all(mode <= xl)) {
      return(xl)

    } else if (all(mode >= xr)) {
      return(xr)

    } else {
      if (f(xl) > f(xl + delta)) {
        new_modes <- c(new_modes, xl)
      }
      if (f(xr) > f(xr - delta)) {
        new_modes <- c(new_modes, xr)
      }
      return(new_modes)

    }
  } else {
    return(new_modes)
  }

}

#' @noRd
check_grid_opt_criteria <- function(symmetric, cnum, dendata) {
  if (dendata$scalable) {
    std_symmetric <- !is.null(symmetric)

    if (cnum %% 2 == 1) {
      scaled_params <- cached_grid_info(cnum + 1)
      which_grid <- "secondary"
    } else {
      scaled_params <- cached_grid_info(cnum)
      which_grid <- "standerd"
    }
    if (!is.null(scaled_params)) {
      if ((std_symmetric &&
           !scaled_params[1]) ||
          (!std_symmetric && scaled_params[1])) {
        msg <- message("you need to delete the ",
                   which_grid,
                   " built-in grid, that has ")
        name <- names(dendata$std_params)
        for (i in (1:(length(scaled_params) - 1))) {
          msg <- message(msg, name[i], " = ", scaled_params[i + 1])

          if (i != length(scaled_params) - 1)
            msg <- message(msg, ", ")
          else
            msg <- message(msg, ".\n")

        }
        stop(msg)
      }
    }
  }
}

#' @noRd
truncate_error_checking <- function(xl, xr, density) {
  stopifnot(
    "xl must be a scaler" = (is.numeric(xl) && length(xl) == 1),
    "xr must be a scaler" = (is.numeric(xr) && length(xr) == 1),
    "xl must be smaller that xr" = (xl < xr),
    "xl must be greater than the density lower bound" = (xl >=  density$lb),
    "xr must be smaller than the density upper bound" = (xr <= density$rb)
  )

  return(list(xl = xl, xr = xr))

}

#' @noRd
grid_error_checking_and_preparation <- function(gp) {
  modes <- gp$target$modes
  f <- gp$target$density
  between_minima <- gp$target$between_minima
  steps <- gp$proposal$steps
  theta <- gp$proposal$pre_acceptance_threshold
  grid_range <- gp$proposal$grid_range
  lb <- gp$target$left_bound
  rb <- gp$target$right_bound

  if (is.null(modes)) {
    stop("Error: 'modes' density modes must be provided.")
  }

  if (!is.function(f)) {
    stop("Error: 'f' density function must be provided.")
  }

  if (!is.null(steps) && steps < 1) {
    stop("Error: 'steps' must be greater than or equal to 1.")
  }

  if (is.null(theta))
    theta <- 0

  if (theta != 0 && !is.null(grid_range)) {
    stop(
      "Error: You must provide either a pre-acceptance threshold 'theta' value or a proposal x-axis limit 'grid_range'."
    )
  }

  if ((theta != 0 || !is.null(grid_range)) && !is.null(steps)) {
    warning(
      "Warning: The pre-acceptance threshold 'theta' value and proposal x-axis limit 'grid_range' will not take effect because you are specifying a target 'steps' number."
    )
  }

  if (theta != 0 && (theta < 0 || theta > 1)) {
    stopifnot(theta >= 0 && theta <= 1,
              "Error: 'theta' must be in the range [0,1]")
  }

  if (!is.null(grid_range)) {
    if (length(grid_range) != 2)
      stop("Error: 'grid_range' must be a vector of two elements.")

    if (grid_range[1] < lb && grid_range[2] > rb)
      stop("Error: 'grid_range' must be within the range of distribution bounds.")

    if (grid_range[1] > modes[1] ||
        grid_range[2] < modes[length(modes)])
      stop("Error: 'grid_range' range must include distribution's modes.")
  } else {
    grid_range <- gp$grid_range <- c(lb, rb)
  }

  if (!is.null(between_minima)) {
    if (between_minima < lb || between_minima > rb)
      stop("Error: 'between_minima' must be within the range of distribution bounds.")

    if (length(between_minima) != (length(modes) - 1))
      stop("Error: 'between_minima' number of minima oints must be equal to number of modes - 1")

    minima_len <- length(between_minima)

    for (i in (1:minima_len)) {
      if (!(between_minima[i] > modes[i] &&
            between_minima[i] < modes[i + 1]))
        stop("Error: 'between_minima' must be inbetween the distribution's modes.")
    }
  }

  gp$proposal$grid_range <- grid_range
  gp$proposal$pre_acceptance_threshold <- theta

  return(gp)

}

#' @import digest digest
#' @noRd
is_valid_grid <- function(grid) {
  if ("lock" %in% names(grid)) {
    temp <- grid[setdiff(names(grid), "lock")]
    key <- digest(temp)

    if (key == grid$lock) {
      return(TRUE)

    }

  }

  return(FALSE)
}


#' @noRd
grid_check_symmetric <- function(gp) {


  if (!is.null(gp$target$symmetric)) {
    # modes <- gp$target$modes
    # rb <- gp$target$right_bound
    # lb <- gp$target$left_bound
    # grid_range <- gp$proposal$grid_range
    # f <- gp$target$density
    # center <- gp$target$symmetric
    # n <- 21
    #
    # if (is.finite(lb) || is.finite(rb)) {
    #   sub <- min(lb, rb)
    # } else {
    #   sub <- 5
    # }
    #
    # vals <- seq(from = center + sub,
    #             to = center - sub,
    #             length.out = n)
    #
    # for (i in as.integer(n / 2)) {
    #   if (f(vals[i]) != f(vals[n - i + 1]))
    #     stop(paste0("the target density is not symmetric around ", center))
    # }
    #
    # gp$target$modes <- modes[modes > center]
    # if (length(modes) == 1) {
    #   modes <- center
    #   grid_range <- c(center, grid_range[2])
    # } else {
    #   if (center %in% modes) {
    #     modes <- modes[modes >= center]
    #   } else {
    #     modes <- modes[modes > center]
    #   }
    # }
    #
    # gp$target$modes <- modes
    # gp$target$left_bound <- center
    # gp$proposal$grid_range <- grid_range
    # gp$target$modes_count <- length(gp$target$modes)
    grid_range <- gp$proposal$grid_range
    center <- gp$target$modes

    gp$target$left_bound <- center
    gp$proposal$grid_range <- c(center, grid_range[2])

  }

  return(gp)

}



#' @noRd
get_buildin_sampling_function <- function(cnum, name) {
  if (cnum %% 2 == 0) {
    even <- TRUE
    cnum_search <- cnum - 1
  } else {
    even <- FALSE
    cnum_search <- cnum
  }



  if (pbgrids[[name]]$c_num == cnum_search) {
    if (even) {
      fun <- function(n) {
        args <- list(n = n)
        return(do.call(paste0(name, "_custom"), args))
      }

    } else {
      fun <- function(n) {
        args <- as.list(c(n = n, pbgrids[[name]]$std_params))
        do.call(paste0(name), args)
      }
    }
  } else {
    stop("Check the grid caching number c_num !")
  }

  return(fun)


}


#' @noRd
built_in_pars_trans <- function(f_params, cnum) {
  for (name in names(pbgrids)) {
    if (pbgrids[[name]]$c_num == cnum ||
        pbgrids[[name]]$c_num + 1 == cnum) {
      return(pbgrids[[name]]$transform_params(f_params))

    }
  }
  return(f_params)
}

#' @noRd
cache_grid_c <- function(c_num, grid) {
  n_params <- length(grid$f_params)

  if (n_params == 0) {
    f_params <- 0
  } else {
    f_params <- built_in_pars_trans(grid$f_params, c_num)

    f_params <- unlist(grid$f_params)
  }


  .Call(
    C_cache_grid,
    c_num,
    grid$grid_data$x,
    grid$grid_data$s_upper,
    grid$grid_data$p_a,
    grid$grid_data$s_upper_lower,
    grid$areas,
    grid$steps_number,
    grid$sampling_probabilities,
    grid$unif_scaler,
    grid$lt_properties,
    grid$rt_properties,
    grid$alpha,
    grid$symmetric,
    f_params,
    n_params,
    grid$grid_bounds[1],
    grid$grid_bounds[2]
  )

}



#' @noRd
cache_user_grid_c <- function(grid) {
  if (!is_valid_grid(grid))
    stop("This grid is not valid")

  if (grid$lock %in% stors_env$user_session_cached_grid_locks[["lock"]]) {
    message("cashed grid already exist, no more cashing cashing!. Just returning c_num !")

    c_num <- stors_env$user_session_cached_grid_locks[stors_env$user_session_cached_grid_locks$lock == grid$lock, ]$cnum

    return(c_num)

  }

  message("this grid has not been cashed before !")

  c_num <- stors_env$user_cnum_counter

  cache_grid_c(c_num, grid)

  user_session_cached_grid_locks <- data.frame(lock = grid$lock, cnum = c_num)

  stors_env$user_session_cached_grid_locks <- rbind(stors_env$user_session_cached_grid_locks,
                                                    user_session_cached_grid_locks)

  stors_env$user_cnum_counter <- stors_env$user_cnum_counter + 1

  return(c_num)
}

#' @noRd
free_cache_cnum_c <- function(c_num) {
  .Call(C_free_cache_cnum, c_num)
}

#' @noRd
save_builtin_grid <- function(c_num, grid) {
  grids_file_path <- file.path(stors_env$builtin_grids_dir, paste0(c_num, ".rds"))
  saveRDS(grid, grids_file_path)
}

#' @noRd
cached_grid_info <- function(cnum) {
  .Call(C_grid_info, cnum)
}


#' Delete Built-in Grids
#'
#' @description
#' This function deletes built-in grids from disk by specifying the sampling function and grid type.
#' It is useful for managing cached grids and freeing up storage space.
#'
#' @param sampling_function String. The name of the sampling distribution's function in STORS.
#' For example, \code{"srgamma"} or \code{"srchisq"}.
#' @param grid_type String. Either \code{"custom"} to delete the custom grid or \code{"scaled"} to delete the scaled grid.
#' Defaults to \code{"custom"}.
#'
#' @details
#' The function looks for the specified grid type associated with the sampling function in the built-in grids directory.
#' If the grid exists, it deletes the corresponding grid file from disk and frees its cached resources.
#' If the specified sampling function or grid type does not exist, an error is thrown.
#'
#' @return
#' A message indicating the status of the deletion process, or an error if the operation fails.
#'
#' @examples
#' # Delete a custom grid for the srgamma function
#' delete_build_in_grid(sampling_function = "srgamma", grid_type = "custom")
#'
#' # Delete a scaled grid for the srnorm function
#' delete_build_in_grid(sampling_function = "srnorm", grid_type = "scaled")
#'
#' @export
delete_build_in_grid <- function(sampling_function, grid_type = "custom") {
  grid_type <- match.arg(grid_type, c("scaled", "custom"))

  if (sampling_function %in% names(pbgrids)) {
    if (grid_type == "scaled") {
      grid_number <- pbgrids[[sampling_function]]$c_num
    } else {
      grid_number <- pbgrids[[sampling_function]]$c_num + 1
    }

  } else {
    stop(paste0("sampling function ", sampling_function, " does not exist!"))

  }

  builtin_grids <- list.files(stors_env$builtin_grids_dir)

  for (grid_name in builtin_grids) {
    grid_path <- file.path(stors_env$builtin_grids_dir, grid_name)
    grid <- readRDS(grid_path)

    if ("lock" %in% names(grid)) {
      temp <- grid[setdiff(names(grid), "lock")]
      key <- digest(temp)

      if (key == grid$lock) {
        if (grid$cnum == grid_number) {
          file.remove(grid_path)
          free_cache_cnum_c(grid$cnum)
          message(" grid number ", grid$cnum, " DELETED !")
        }

      }

    }

  }

}

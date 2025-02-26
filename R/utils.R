#' @noRd
create_function <- function(fun, arg_list) {

  fun_args <- names(formals(fun))

  arg_list <- arg_list[names(arg_list) %in% fun_args]

  build_closure <- function(fun, ...) {
    return(function(x) {
      fun(x, ...)
    })
  }

  f <- do.call("build_closure", c(fun, arg_list))

  return(f)
}



#' @noRd
fix_function <- function(func) {
  func_env <- environment(func)
  func_body <- body(func)
  func_args <- names(formals(func))

  if (length(func_args) != 1) stop("all provided function must has only one parameter")

  # Remove parameters from the environment to exclude them from substitution
  captured_env <- as.list(func_env)
  captured_env <- captured_env[!names(captured_env) %in% func_args]

  # Replace free variables in the body with their values from the environment
  fixed_body <- eval(substitute(substitute(func_body, captured_env), list(func_body = func_body)))

  # Create a new function with the modified body
  new_func <- func
  body(new_func) <- fixed_body
  environment(new_func) <- baseenv()  # Reset the environment to avoid referencing external variables

  new_func
}


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
check_proposal_opt_criteria <- function(symmetric, cnum, dendata) {
  if (dendata$scalable) {
    std_symmetric <- !is.null(symmetric)

    if (cnum %% 2 == 1) {
      scaled_params <- cached_proposal_info(cnum + 1)
      which_proposal <- "secondary"
    } else {
      scaled_params <- cached_proposal_info(cnum)
      which_proposal <- "standerd"
    }
    if (!is.null(scaled_params)) {
      if ((std_symmetric &&
           !scaled_params[1]) ||
          (!std_symmetric && scaled_params[1])) {
        msg <- message("you need to delete the ",
                   which_proposal,
                   " built-in proposal, that has ")
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
    "xl must be greater than the density lower bound" = (xl >=  density$lower),
    "xr must be smaller than the density upper bound" = (xr <= density$upper)
  )

  return(list(xl = xl, xr = xr))

}

#' @noRd
proposal_error_checking_and_preparation <- function(gp) {
  modes <- gp$target$modes
  f <- gp$target$density
  between_minima <- gp$target$between_minima
  steps <- gp$proposal$steps
  theta <- gp$proposal$pre_acceptance_threshold
  proposal_range <- gp$proposal$proposal_range
  lower <- gp$target$left_bound
  upper <- gp$target$right_bound

  if (is.null(modes)) {
    stop("Error: 'modes' density modes must be provided.")
  }

  if (!is.function(f)) {
    stop("Error: 'f' density function must be provided.")
  }

  if (!is.null(steps) && steps < 1) {
    stop("Error: 'steps' must be greater than or equal to 1.")
  }

  # if (theta != 0 && !is.null(proposal_range)) {
  #   stop(
  #     "Error: You must provide either a pre-acceptance threshold 'theta' value or a proposal x-axis limit 'proposal_range'."
  #   )
  # }
  #
  # if ((theta != 0 || !is.null(proposal_range)) && !is.null(steps)) {
  #   warning(
  #     "Warning: The pre-acceptance threshold 'theta' value and proposal x-axis limit 'proposal_range' will not take effect because you are specifying a target 'steps' number."
  #   )
  # }

  if ( !is.numeric(theta) || (theta <  0.1 || theta > 1)) {
    stopifnot(theta >= 0 && theta <= 1,
              "Error: 'theta' must be in the range [0.1,1]")
  }

  if (!is.null(proposal_range)) {
    if (length(proposal_range) != 2)
      stop("Error: 'proposal_range' must be a vector of two elements.")

    if (proposal_range[1] < lower && proposal_range[2] > upper)
      stop("Error: 'proposal_range' must be within the range of distribution bounds.")

    if (proposal_range[1] > modes[1] ||
        proposal_range[2] < modes[length(modes)])
      stop("Error: 'proposal_range' range must include distribution's modes.")
  } else {
    proposal_range <- gp$proposal_range <- c(lower, upper)
  }

  if (!is.null(between_minima)) {
    if (between_minima < lower || between_minima > upper)
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

  gp$proposal$proposal_range <- proposal_range
  gp$proposal$pre_acceptance_threshold <- theta

  return(gp)

}

#' @import digest digest
#' @noRd
is_valid_proposal <- function(proposal) {
  if ("lock" %in% names(proposal)) {
    temp <- proposal[setdiff(names(proposal), "lock")]
    key <- digest(temp)

    if (key == proposal$lock) {
      return(TRUE)

    }

  }

  return(FALSE)
}


#' @noRd
proposal_check_symmetric <- function(gp) {


  if (!is.null(gp$target$symmetric)) {
    # modes <- gp$target$modes
    # upper <- gp$target$right_bound
    # lower <- gp$target$left_bound
    # proposal_range <- gp$proposal$proposal_range
    # f <- gp$target$density
    # center <- gp$target$symmetric
    # n <- 21
    #
    # if (is.finite(lower) || is.finite(upper)) {
    #   sub <- min(lower, upper)
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
    #   proposal_range <- c(center, proposal_range[2])
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
    # gp$proposal$proposal_range <- proposal_range
    # gp$target$modes_count <- length(gp$target$modes)
    proposal_range <- gp$proposal$proposal_range
    center <- gp$target$modes

    gp$target$left_bound <- center
    gp$proposal$proposal_range <- c(center, proposal_range[2])

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



  if (built_in_proposals[[name]]$c_num == cnum_search) {
    if (even) {
      fun <- function(n) {
        args <- list(n = n)
        return(do.call(paste0(name, "_custom"), args))
      }

    } else {
      fun <- function(n) {
        args <- as.list(c(n = n, built_in_proposals[[name]]$std_params))
        do.call(paste0(name), args)
      }
    }
  } else {
    stop("Check the proposal caching number c_num !")
  }

  return(fun)


}


#' @noRd
built_in_pars_trans <- function(f_params, cnum) {
  for (name in names(built_in_proposals)) {
    if (built_in_proposals[[name]]$c_num == cnum ||
        built_in_proposals[[name]]$c_num + 1 == cnum) {
      return(built_in_proposals[[name]]$transform_params(f_params))

    }
  }
  return(f_params)
}

#' @noRd
cache_proposal_c <- function(c_num, proposal) {
  n_params <- length(proposal$f_params)

  f_params <- 0

  if (n_params > 0) {
    f_params <- built_in_pars_trans(proposal$f_params, c_num)

    f_params <- unlist(proposal$f_params)
  }

  .Call(
    C_cache_grid,
    c_num,
    proposal$data$x,
    proposal$data$s_upper,
    proposal$data$p_a,
    proposal$data$s_upper_lower,
    proposal$areas,
    proposal$steps_number,
    proposal$sampling_probabilities,
    proposal$unif_scaler,
    proposal$lt_properties,
    proposal$rt_properties,
    proposal$alpha,
    proposal$symmetric,
    f_params,
    n_params,
    proposal$proposal_bounds[1],
    proposal$proposal_bounds[2]
  )

}



#' @noRd
cache_user_proposal_c <- function(proposal) {
  if (!is_valid_proposal(proposal))
    stop("This proposal is not valid")

  if (proposal$lock %in% stors_env$user_session_cached_proposals_locks[["lock"]]) {
     message("Detected cached proposal for this distribution already ... replacing with new proposal.")

    c_num <- stors_env$user_session_cached_proposals_locks[stors_env$user_session_cached_proposals_locks$lock == proposal$lock, ]$cnum

    return(c_num)

  }

  if(stors_env$user_cnum_counter > 199){
    stop("Cache limit exceeded: Users can cache up to 100 proposals in memory.")
  }

  c_num <- stors_env$user_cnum_counter

  cache_proposal_c(c_num, proposal)

  user_session_cached_proposals_locks <- data.frame(lock = proposal$lock, cnum = c_num)

  stors_env$user_session_cached_proposals_locks <- rbind(stors_env$user_session_cached_proposals_locks,
                                                    user_session_cached_proposals_locks)

  stors_env$user_cnum_counter <- stors_env$user_cnum_counter + 1

  return(c_num)
}

#' @noRd
free_cache_cnum_c <- function(c_num) {
  .Call(C_free_cache_cnum, c_num)
}

#' @noRd
save_builtin_proposal <- function(c_num, proposal) {
  proposals_file_path <- file.path(stors_env$builtin_proposals_dir, paste0(c_num, ".rds"))
  saveRDS(proposal, proposals_file_path)
}

#' @noRd
cached_proposal_info <- function(cnum) {
  .Call(C_grid_info, cnum)
}


#' Delete Built-in Proposal
#'
#' @description
#' This function deletes built-in proposals from disk by specifying the sampling function and proposal type.
#' It is useful for managing cached proposals and freeing up storage space.
#'
#' @param sampling_function String. The name of the sampling distribution's function in STORS.
#' For example, \code{"srgamma"} or \code{"srchisq"}.
#' @param proposal_type String. Either \code{"custom"} to delete the custom proposal or \code{"scaled"} to delete the scaled proposal.
#' Defaults to \code{"custom"}.
#'
#' @details
#' The function looks for the specified proposal type associated with the sampling function in the built-in proposals directory.
#' If the proposal exists, it deletes the corresponding proposal file from disk and frees its cached resources.
#' If the specified sampling function or proposal type does not exist, an error is thrown.
#'
#' @return
#' A message indicating the status of the deletion process, or an error if the operation fails.
#'
#' @examples
#'
#' # The following examples are commented, since if they are run the srgamma and
#' # srnorm samplers will no longer work until a new grid is built for them.
#' # This causes problems if the examples are run by CRAN checks or the website
#' # build system, hence they are commented.
#'
#' # Delete a custom proposal for the srgamma function (uncomment to run)
#' #delete_build_in_proposal(sampling_function = "srgamma", proposal_type = "custom")
#'
#' # Delete a scaled proposal for the srnorm function (uncomment to run)
#' #delete_build_in_proposal(sampling_function = "srnorm", proposal_type = "scaled")
#'
#' @export
delete_build_in_proposal <- function(sampling_function, proposal_type = "custom") {
  proposal_type <- match.arg(proposal_type, c("scaled", "custom"))

  if (sampling_function %in% names(built_in_proposals)) {
    if (proposal_type == "scaled") {
      proposal_number <- built_in_proposals[[sampling_function]]$c_num
    } else {
      proposal_number <- built_in_proposals[[sampling_function]]$c_num + 1
    }

  } else {
    stop(paste0("sampling function ", sampling_function, " does not exist!"))

  }

  builtin_proposals <- list.files(stors_env$builtin_proposals_dir)

  for (proposal_name in builtin_proposals) {
    proposal_path <- file.path(stors_env$builtin_proposals_dir, proposal_name)
    proposal <- readRDS(proposal_path)

    if ("lock" %in% names(proposal)) {
      temp <- proposal[setdiff(names(proposal), "lock")]
      key <- digest(temp)

      if (key == proposal$lock) {
        if (proposal$cnum == proposal_number) {
          file.remove(proposal_path)
          free_cache_cnum_c(proposal$cnum)
          message(" proposal number ", proposal$cnum, " DELETED !")
        }

      }

    }

  }

}

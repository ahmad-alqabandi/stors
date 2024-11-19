#' @noRd
load_ddls_name = function(dist_name, is_symmetric = FALSE) {
  if (is_symmetric) {
    assign("C_sampling_fun", get(paste0("C_", dist_name, "_sym")), envir = environment(get(dist_name)))
  } else{
    assign("C_sampling_fun", get(paste0("C_", dist_name)), envir = environment(get(dist_name)))
  }
  
}

#' @noRd
check_grid_optimization_criteria = function(symmetric, cnum, dendata) {
  if (dendata$scalable) {
    std_symmetric = !is.null(symmetric)
    
    if (cnum %% 2 == 1) {
      scaled_params <- cached_grid_info(cnum + 1)
      which_grid <- "secondary"
    } else{
      scaled_params <- cached_grid_info(cnum)
      which_grid <- "standerd"
    }
    if (!is.null(scaled_params)) {
      if ((std_symmetric &&
           !scaled_params[1]) || (!std_symmetric && scaled_params[1])) {
        msg <- cat("you need to delete the ",
                   which_grid,
                   " built-in grid, that has ")
        name <- names(dendata$std_params)
        for (i in (1:(length(scaled_params) - 1)))
        {
          msg <- cat(msg, name[i], " = ", scaled_params[i + 1])
          if (i != length(scaled_params) - 1)
            msg <- cat(msg, ", ")
          else
            msg <- cat(msg, ".\n")
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
grid_error_checking_and_preparation = function(gp) {
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
  } else{
    grid_range <- gp$grid_range <- c(lb, rb)
  }
  
  if (!is.null(between_minima)) {
    if (between_minima < lb || between_minima > rb)
      stop("Error: 'between_minima' must be within the range of distribution bounds.")
    
    if (length(between_minima) != (length(modes) - 1))
      stop("Error: 'between_minima' number of minima oints must be equal to number of modes - 1")
    
    minima_len = length(between_minima)
    
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
is_valid_grid = function(grid) {
  stopifnot(
    " This grid is not optmized using build_grid() " = digest(grid) %in% stors_env$created_girds_Id
  )
  
}

#' @noRd
grid_check_symmetric <- function(gp) {
  if (!is.null(gp$target$symmetric)) {
    # if( is.null(gp$target$symmetric_around_value) ){
    #   if(gp$target$modes_count > 1){
    #     stop("you need to provide symmetric_around_value to build symmetric grid.")
    #   }
    #   warning("symmetric_around_value has been set equale to the distrebution's mode.")
    #   gp$target$symmetric_around_value <- gp$target$modes
    # }
    
    modes <- gp$target$modes
    rb <- gp$target$right_bound
    lb <- gp$target$left_bound
    grid_range <- gp$proposal$grid_range
    
    f <- gp$target$density
    
    center <- gp$target$symmetric
    
    n <- 21
    
    if (is.finite(lb) || is.finite(rb)) {
      sub <- min(lb, rb)
    } else{
      sub <- 5
    }
    
    vals <- seq(from = center + sub,
                to = center - sub ,
                length.out = n)
    
    for (i in as.integer(n / 2)) {
      if (f(vals[i]) != f(vals[n - i + 1]))
        stop(paste0("the target density is not symmetric around ", center))
    }
    
    gp$target$modes <- modes[modes > center]
    if (length(modes) == 1) {
      modes <- center
      grid_range <- c(center, grid_range[2])
    } else{
      if (center %in% modes) {
        modes <- modes[modes >= center]
      } else{
        modes <- modes[modes > center]
      }
    }
    
    gp$target$modes <- modes
    gp$target$left_bound <- center
    gp$proposal$grid_range <- grid_range
    gp$target$modes_count <- length(gp$target$modes)
  }
  
  return(gp)
  
}


#' @noRd
cache_grid_c <- function(Cnum, grid) {
  n_params <- length(grid$f_params)
  
  if (n_params == 0) {
    f_params <- 0
  } else{
    f_params <- unlist(grid$f_params)
  }
  
  .Call(
    C_cache_grid,
    Cnum,
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
    n_params
  )
  
}

#' @noRd
free_cache_cnum_c <- function(Cnum) {
  .Call(C_free_cache_cnum, Cnum)
}

#' @noRd
save_builtin_grid <- function(Cnum, grid) {
  grids_file_path <- file.path(stors_env$user_dirs$builtin_dir, paste0(Cnum, ".rds"))
  
  saveRDS(grid, grids_file_path)
}

#' @noRd
delete_build_in_grid_cnum <- function(Cnum) {
  path <- file.path(stors_env$user_dirs$builtin_dir, paste0(Cnum, ".rds"))
  if (file.exists(path)) {
    file.remove(path)
    
    
    
  } else{
    stop("This grid does not exist.\n")
  }
  
  free_cache_cnum_c(Cnum)
  
}

#' @export
cached_grid_info = function(cnum) {
  .Call(C_grid_info, cnum)
}

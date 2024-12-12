
#' @export
srbeta_custom <- function(n = 1, x = NULL) {
  .Call(C_srbeta_custom_check, n, x)
}





#' @export
srbeta_optimize = function(
    shape1 = NULL,
    shape2 = NULL,
    xl = 0,
    xr = 1,
    steps = NULL,
    grid_range = NULL,
    theta = NULL,
    target_sample_size = 1000,
    verbose = FALSE) {
  
  symmetric <- NULL
  
  dist_name <- 'srbeta'
  
  dendata <- pbgrids[[dist_name]]
  
  if (shape1 <= 1 || shape2 <= 1) {
    message("Grid building is not available for shape1 <= 1 or shape2 <= 1.")
    return()
  }
  
  f_params <- list(shape1 = shape2, shape2 = shape2)
  
  if(dendata$scalable){
    
    isnull <- sapply(f_params, is.null)
    
    if(all(isnull)){
      cnum <- dendata$Cnum
      grid_type = "scaled"
    }else{
      cnum <- dendata$Cnum + 1
      grid_type = "custom"
    }
    
    f_params <- ifelse(isnull, dendata$std_params, f_params)
    
  }else{
    cnum <- dendata$Cnum + 1
    grid_type = "custom"
  }
  
  modes <- dendata$set_modes(f_params$shape1, f_params$shape2)
  
  f <- dendata$create_f(f_params$shape1, f_params$shape2)
  
  check_grid_optimization_criteria(symmetric, cnum, dendata)
  
  grid_optimizer(dendata, dist_name, xl, xr, f, modes, f_params, steps,
                 grid_range, theta, target_sample_size,
                 grid_type, symmetric, cnum, verbose)
  
}
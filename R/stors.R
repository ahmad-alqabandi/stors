dyn.load("./src/stors.so")


#' stros for sampling
#'
#' @param n sample size
#' @param grid proposal grid
#' @param f target density
#' @return sample of size n
#' @export
#'
stors <- function( dist_name) {

  data_dir <- tools::R_user_dir("stors", "data")
  
  stopifnot(`there are no grids in user_data directory` = dir.exists(data_dir))
  
  grid_dir = file.path(data_dir, paste0(dist_name,".rds"))
  
  stopifnot(`there is no grid for 'grid_dir' density` = file.exists(grid_dir))
  
  grid = readRDS(grid_dir)
  
  grid$dens_func = eval(parse(text=grid$dens_func))
  
  function(n){
    
  .Call("stors" ,
        n,
        grid$grid_data$x,
        grid$grid_data$s_upper_lower,
        grid$grid_data$p_a,
        grid$steps_number,
        grid$sampling_probabilities,
        grid$unif_scaler,
        grid$grid_data$s_upper,
        grid$lt_properties,
        grid$rt_properties,
        grid$dens_func,
        new.env(),
        PACKAGE = "stors")
    
  }
}



#' Title
#'
#' @param n 
#' @param grid 
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
stors2 <- function( n, grid, f) {
  
    .Call("stors" ,
          n,
          grid$grid_data$x,
          grid$grid_data$s_upper_lower,
          grid$grid_data$p_a,
          grid$steps_number,
          grid$sampling_probabilities,
          grid$unif_scaler,
          grid$grid_data$s_upper,
          grid$lt_properties,
          grid$rt_properties,
          f,
          new.env(),
          PACKAGE = "stors")
    
  
}


#' Title
#'
#' @param n 
#' @param x 
#' @param lower 
#' @param pa 
#' @param steps 
#' @param pro 
#' @param unis 
#' @param upper 
#' @param lt 
#' @param rt 
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
stors3 <- function(n, x,lower, pa,steps,pro,unis, upper,lt,rt, f) {
  
  .Call("stors" ,
        n,
       x,
       lower,
        pa,
        steps,
        pro,
       unis,
       upper,
        lt,
        rt,
        f,
        new.env(),
        PACKAGE = "stors")
  
  
}



# 
# grid_size <- function( obj){
#   .Call("print_size", obj, PACKAGE = "stors")
# }
# 
# 
# 
# pre_fetch <- function( obj, size){
#   invisible(.Call("pre_fetch", obj, size, PACKAGE = "stors"))
# }

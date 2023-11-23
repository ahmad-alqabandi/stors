#' @useDynLib stors, .registration = TRUE, .fixes = "C_"

stors_env <- new.env(parent = emptyenv())

.onLoad <- function(lib, pkg) {
  
  data_dir <- tools::R_user_dir("stors", "data")
  
  builtin_dir = file.path(data_dir,"biultin_grids")
  

  if (!dir.exists(builtin_dir)) dir.create(builtin_dir, recursive = TRUE)

  stors_grids_path <- file.path(builtin_dir, "grids.rds")
  
  
  if (file.exists(stors_grids_path)) {
    
    grids <- readRDS(stors_grids_path)

    for (name in grids$biultin$names) {
      
      if (grids$biultin[[name]]$opt) {
        opt_grid <- readRDS(file.path(stors_grids_path, paste0(grids$biultin[[name]]$Cnum, ".rds")))
        cash_grid_c(grids$biultin[[name]]$Cnum, opt_grid)
      }
      
    }
    
  } else{
    grids <- list(biultin = list(names = c("srnorm", "sgamma"),
                          srnorm = list(opt = FALSE, Cnum = 0),
                          sgamma = list(opt = FALSE, Cnum = 1),
                          builtin_num = 2),
                  user = data.frame( name = character(), efficiency = double()) )
  } 
  
  user_cached_grids <- data.frame( Id = character(), Cnum = integer())
  user_dirs <- list(data_dir = data_dir, builtin_dir = builtin_dir)
  created_girds_Id = character()
    
  
  assign("grids", grids, envir = stors_env)
  assign("user_cached_grids", user_cached_grids, envir = stors_env)
  assign("user_dirs", user_dirs, envir = stors_env)
  assign("created_girds_Id", created_girds_Id, envir = stors_env)

}



.onUnload <- function(...) {
  
  stors_env_path <- file.path(stors_env$user_dirs$builtin_dir, "grids.rds")

  saveRDS(stors_env$grids, stors_env_path)

  .Call(C_free_cache)
}

#' @useDynLib stors, .registration = TRUE, .fixes = "C_"

# pbgrids : 1 Normal

pbgrids <- list(
  srnorm = list(
    Cnum = 1,
    tails_method = "ARS",
    f = function(x) {
      0.3989423 * exp(-0.5 * x * x)
    },
    h = function(x) {
      log(0.3989423) - (x * x) * (1 / 2)
    },
    h_prime = function(x) {
      -x
    },
    modes = 0,
    lb = -Inf,
    rb = Inf,
    is.lb = -8.125891,
    is.rb = 8.125891
  ),
  srlaplace = list(
    Cnum = 2,
    tails_method = "IT",
    f = function(x) {
      0.5 * exp(-abs(x))
    },
    cdf = function(x) {
      if (x <= 0) {
        0.5 * exp(x)
      } else {
        1 - 0.5 * exp(-x)
      }
    },
    modes = 0,
    lb = -Inf,
    rb = Inf,
    is.lb = -16,
    is.rb = 16
  ),
  srexp = list(
    Cnum = 3,
    tails_method = "IT",
    f = function(x) {
      exp(-x)
    },
    cdf = function(x) {
      1 - exp(-x)
    },
    modes = 0,
    lb = 0,
    rb = Inf,
    is.lb = 0,
    is.rb = 36.04365
  )
  
)




stors_env <- new.env(parent = emptyenv())

.onLoad <- function(lib, pkg) {
  # look at all exported functions from STORS namespace
  # get all functions start with "sr"
  # add them as string in builtin list
  
  data_dir <- tools::R_user_dir("stors", "data")
  
  builtin_dir = file.path(data_dir, "builtin_grids")
  
  
  if (!dir.exists(builtin_dir))
    dir.create(builtin_dir, recursive = TRUE)
  
  
  stors_grids_path <- file.path(builtin_dir, "grids.rds")
  
  
  if (file.exists(stors_grids_path)) {
    grids <- readRDS(stors_grids_path)
    
    for (name in grids$builtin$names) {
      if (grids$builtin[[name]]$opt) {
        opt_grid <- readRDS(file.path(builtin_dir, paste0(grids$builtin[[name]]$Cnum, ".rds")))
        cache_grid_c(grids$builtin[[name]]$Cnum, opt_grid)
      }
      
    }
    
    first_time_load = FALSE
    
  } else{
    grids <- list(
      builtin = list(
        names = names(pbgrids),
        builtin_num = length(pbgrids)
      ),
      user = data.frame(name = character(), efficiency = double())
    )
    
    
    for (name in names(pbgrids)) {
      grids$builtin[[name]] = list(opt = FALSE, Cnum = pbgrids[[name]]$Cnum)
    }
    
    first_time_load = TRUE
    
  }
  
  user_cached_grids <- data.frame(Id = character(), Cnum = integer())
  user_dirs <- list(data_dir = data_dir, builtin_dir = builtin_dir)
  created_girds_Id = character()
  
  assign("grids", grids, envir = stors_env)
  assign("user_cached_grids", user_cached_grids, envir = stors_env)
  assign("user_dirs", user_dirs, envir = stors_env)
  assign("created_girds_Id", created_girds_Id, envir = stors_env)
  
  
  if (first_time_load) {
    for (name in names(pbgrids)) {
      grid_optimizer(density_name = name, steps = 4091)
    }
  }
  
}



.onUnload <- function(...) {
  stors_env_path <- file.path(stors_env$user_dirs$builtin_dir, "grids.rds")
  
  saveRDS(stors_env$grids, stors_env_path)
  
  .Call(C_free_cache)
}

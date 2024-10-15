#' @useDynLib stors, .registration = TRUE, .fixes = "C_"


pbgrids <- list(
  srnorm = list(
    Cnum = 1,
    tails_method = "ARS",
    symmetric = FALSE,
    create_f = function(mu, sd) {
      function(x)
        ((1.0 / (sd * 2.50662827463) *
            exp(-0.5 * ((
              x - mu
            ) / sd) *
              ((
                x - mu
              ) / sd))))
    },
    set_modes = function(mu = 0)
      mu,
    lb = -Inf,
    rb = Inf
  ),
  srlaplace = list(
    Cnum = 2,
    symmetric = FALSE,
    tails_method = "IT",
    create_f = function(mu, b) {
      function(x) {
        (1 / (2 * b)) * exp(-abs(x - mu) / b)
      }
    },
    create_cdf = function(mu, b) {
      function(x) {
        if (x <= mu) {
          0.5 * exp((x - mu) / b)
        } else {
          1 - 0.5 * exp(-(x - mu) / b)
        }
      }
    },
    set_modes = function(mu = 0)
      mu,
    lb = -Inf,
    rb = Inf
  ),
  old_srnorm = list(
    Cnum = 3,
    tails_method = "ARS",
    symmetric = FALSE,
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
    rb = Inf
  ),
  srnorm_symmetric = list(
    Cnum = 4,
    tails_method = "ARS",
    symmetric = TRUE,
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
    rb = Inf
  ),
  srchisq = list(
    Cnum = 5,
    tails_method = "ARS",
    symmetric = FALSE,
    create_f = function(df) {
      function(x)
      {
        fx <- (1 / (2 ^ (df / 2) * gamma(df / 2))) * x ^ (df / 2 - 1) * exp(-x / 2)
        if(is.nan(fx)) return(0) else return(fx) # NOTE: this precision issue is not solved for in the C code
      }
    },
    set_modes = function(df = 1) {
      max(df - 2, 0)
    },
    lb = 0,
    rb = Inf
  ), srgamma = list(
    Cnum = 6,
    tails_method = "ARS",
    symmetric = FALSE,
    create_f = function(shape = 1, rate = 1, scale = 1/rate) {
      function(x)
      {
        if (x < 0) {
          return(0)
        }
        fx <- (1 / (gamma(shape) * scale ^ shape) * x ^ (shape - 1) * exp(-x / scale))
        if(is.nan(fx)) return(0) else return(fx) # NOTE: this precision issue is not solved for in the C code
        
      }
    },
    set_modes = function(shape, scale) {
      if(shape < 1) 0 else (shape-1) * scale
    },
    lb = 0,
    rb = Inf
  )
)

stors_env <- new.env(parent = emptyenv())

.onLoad <- function(lib, pkg) {
  
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
      fun_text <- paste0(name,'_optimize()')
      fun_parse <- parse(text = fun_text)
      eval(fun_parse)
    }
  }
  
}



.onUnload <- function(...) {
  stors_env_path <- file.path(stors_env$user_dirs$builtin_dir, "grids.rds")
  
  saveRDS(stors_env$grids, stors_env_path)
  
  .Call(C_free_cache)
}

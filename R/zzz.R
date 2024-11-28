#' @useDynLib stors, .registration = TRUE, .fixes = "C_"

# Globals:


pbgrids <- list(
  srnorm = list(
    Cnum = 1,
    tails_method = "ARS",
    scalable = TRUE,
    std_params = list(mean = 0, sd = 1),
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
    Cnum = 3,
    tails_method = "IT",
    scalable = TRUE,
    std_params = c(0,1),
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
  srexp = list(
    Cnum = 5,
    tails_method = "IT",
    scalable = TRUE,
    std_params = c(1),
    create_f = function(rate) {
      function(x)
      {
        fx <- rate * exp(-rate * x)
        fx <- ifelse(is.nan(fx), 0, fx)
        return(fx)
      }
    },
    create_cdf = function(rate) {
      function(x) {
        return(1 - exp(- rate * x))
      }
      },
    lb = 0,
    rb = Inf
  )
  ,
  srchisq = list(
    Cnum = 7,
    tails_method = "ARS",
    scalable = FALSE,
    create_f = function(df) {
      function(x)
      {
        fx <- (1 / (2 ^ (df / 2) * gamma(df / 2))) * x ^ (df / 2 - 1) * exp(-x / 2)
        fx <- ifelse(is.nan(fx), 0, fx)
        return(fx)
      }
    },
    set_modes = function(df = 1) {
      max(df - 2, 0)
    },
    lb = 0,
    rb = Inf
  ), srgamma = list(
    Cnum = 9,
    tails_method = "ARS",
    scalable = FALSE,
    create_f = function(shape = 1, rate = 1, scale = 1/rate) {
      function(x)
      {
        if (x < 0) {
          return(0)
        }
        fx <- (1 / (gamma(shape) * scale ^ shape) * x ^ (shape - 1) * exp(-x / scale))
        fx <- ifelse(is.nan(fx), 0, fx)
        return(fx)
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
  
  builtin_grids_dir <- file.path(data_dir, "builtin_grids")
  
  user_grids_dir <- file.path(data_dir, "user_grids")
  
  
  if (!dir.exists(builtin_grids_dir))
    dir.create(builtin_grids_dir, recursive = TRUE)
  
  
  if (!dir.exists(user_grids_dir))
    dir.create(user_grids_dir, recursive = TRUE)
  
  created_girds_Id = character()
  
  assign("builtin_grids_dir", builtin_grids_dir, envir = stors_env)
  assign("user_grids_dir", user_grids_dir, envir = stors_env)
  assign("created_girds_Id", created_girds_Id, envir = stors_env)
  
  
  #  stors_grids_path <- file.path(builtin_grids_dir, "grids.rds")
  
  builtin_grids <- list.files(builtin_grids_dir)
  
  if(length(builtin_grids) == 0 ){
    
    # here we have to optimize for all scalable grids
    
    # for (name in names(pbgrids)) {
    
    name <- 'srnorm' # temp
    
      if(pbgrids[[name]]$scalable){
        fun_name <- paste0(name,'_optimize')
        do.call(fun_name, args = pbgrids[[name]]$std_params)
        
      }

    # }
    
  }else{
    
    #here we have to load all RDS files , check validation using digest
    # cache all grids
    
    for(grid_name in builtin_grids){
      
      grid_path <- file.path(builtin_grids_dir, grid_name)
      grid <- readRDS(grid_path)
      
      if("lock" %in% names(grid)){
        
        temp <-grid[setdiff(names(grid),"lock")]
        key <- digest(temp)
        
        if( key == grid$lock){
          
          cat(" grid number ", grid$cnum, " CACHED !")
          cache_grid_c(grid$cnum, grid)
          
        }
        
      }
      
    }
    
  }
  
}



.onUnload <- function(...) {
  # stors_env_path <- file.path(stors_env$user_dirs$builtin_grids_dir, "grids.rds")
  
  # saveRDS(stors_env$grids, stors_env_path)
  
  .Call(C_free_cache)
}

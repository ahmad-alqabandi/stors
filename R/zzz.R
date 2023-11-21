#' @useDynLib stors, .registration = TRUE, .fixes = "C_"

grids_env <- new.env(parent = emptyenv())


.onLoad <- function(lib, pkg) {
  data_dir <- tools::R_user_dir("stors", "data")

  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  cache_dir <- tools::R_user_dir("stors", "cache")

  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)


  grids_env_path <- file.path(data_dir, "grids_config.rds")

  if (file.exists(grids_env_path)) {
    grids_config <- readRDS(grids_env_path)

    for (name in grids_config$names) {
      if (grids_config[[name]]$opt) {
        opt_grid <- readRDS(file.path(grids_config$data_dir, paste0(grids_config[[name]]$Cnum, ".rds")))
        cash_grid_c(grids_config[[name]]$Cnum, opt_grid)
      }
    }
  } else {
    grids_config <- list(
      names = c("snorm", "sgamma"),
      snorm = list(opt = FALSE, Cnum = 0),
      sgamma = list(opt = FALSE, Cnum = 1),
      builtin_num = 2,
      grids = data.frame(name = character(), Id = character(), Cnum = integer()),
      data_dir = data_dir,
      cache_dir = cache_dir,
      creatd_Id = character()
    )
  }

  assign("grids_config", grids_config, envir = grids_env)
}



.onUnload <- function(...) {
  grids_env_path <- file.path(grids_env$grids_config$data_dir, "grids_config.rds")

  saveRDS(grids_env$grids_config, grids_env_path)

  .Call(C_free_cache)
}

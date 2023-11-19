#' @useDynLib stors, .registration = TRUE, .fixes = "C_"

grids_env <- new.env(parent = emptyenv())

grids <- data.frame(Id = c(NA, NA), Cnum = c(1, 2))

assign("grids", grids, envir = grids_env)


.onLoad <- function(lib, pkg) {
  # library.dynam("stors", pkg, lib )
}


onUnload <- function(libpath) {
  # if (is.loaded("stors", PACKAGE="stors")) {
  #   library.dynam.unload("stors", libpath)
  # }
}

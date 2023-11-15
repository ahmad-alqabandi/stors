path <- "./R/grid_pre_builds/grids.Rdata"


if (file.exists(path) == FALSE) {


} else {
  load("./R/grid_pre_builds/grids.Rdata")
}


# save R data for testing

grid_builder(lb = -Inf, rb = Inf, mode = 0, f = f_snorm, h = h_snorm, h_prime = h_prime_snorm)
 
 
save(grid, file = "./R/grid_pre_builds/grids.Rdata")


file()

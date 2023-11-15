modes <- 0

f_norm <- function(x) {
  0.3989423 * exp(-0.5 * x * x)
}


h_norm <- function(x) {
  log(0.3989423) - (x * x) * (1 / 2)
}


h_prime_norm <- function(x) {
  -x
}


set_grid(lb = -Inf, rb = Inf, mode = modes, f = f_norm, h = h_norm, h_prime = h_prime_norm, dist_name = "normal", to = "set" )

grids_list_print()

grid = stors::read_grid(dist_name = "normal")

snorm = stors::stors("normal")

x =  grid$grid_data$x
lower = grid$grid_data$s_upper_lower
pa = grid$grid_data$p_a
steps = grid$steps_number
pro = grid$sampling_probabilities
unis= grid$unif_scaler
upper = grid$grid_data$s_upper
lt = grid$lt_properties
rt = grid$rt_properties


n=1000

microbenchmark::microbenchmark(
  snorm(n),
  stors::stors2(n,grid, f_norm),
  stors::stors3(n, x, lower, pa, steps, pro, unis, upper, lt, rt, f_norm),
  rnorm(n),
  times = 10000
)

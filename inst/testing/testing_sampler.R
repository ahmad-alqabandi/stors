
##    STORS


modes = 0

f_norm = function(x) {
  0.3989423 * exp(-0.5 * x * x)
}

h_norm = function(x) {
  log(0.3989423) - (x * x) * (1 / 2)
}

h_prime_norm = function(x) {
  -x
}

norm_grid = build_grid(lb = -Inf, rb = Inf, modes , f = f_norm, h = h_norm, h_prime = h_prime_norm,
                       steps = 1000, verbose = TRUE,  target_sample_size = 10000)

plot(norm_grid)

? save_grid

? print_grids

? delete_grid

? load_grid

norm_grid

save_grid(norm_grid, "normal")

print_grids()

steps = 2000

norm_grid = build_grid(lb = -Inf, rb = Inf, modes , f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = steps)

norm_grid

plot(norm_grid)

norm_grid

norm_sampler = stors(norm_grid)

hist(norm_sampler(1000000))

norm_trunc = stors(norm_grid, -2,2)

hist(norm_trunc(10000))

#======================= Built-in densities

modes_multi = c(0.00134865,3.99865)

p <- 0.5

q <- 0.5

f_multi <- function(x) {
  0.5 * (sqrt(2 * pi))^(-1) * exp(-(x^2)/2) + 0.5 * ( sqrt(2 * pi))^(-1) * exp(-((x - 4)^2)/2)
}

h_multi <- function(x) {
  log(  p * (sqrt(2 * pi))^(-1) * exp(-(x^2)/2) + q * (sqrt(2 * pi))^(-1) * exp(-((x - 4)^2)/2))
}

h_prime_multi <- function(x) {
  (-(exp(-1/2 * (-4 + x)^2) * q * (-4 + x))/sqrt(2 * pi) - (exp(-x^2/2) * p * x)/sqrt(2 * pi))/((exp(-x^2/2) * p)/sqrt(2 * pi) + (exp(-1/2 * (-4 + x)^2) * q)/sqrt(2 * pi))
}

multi_grid = build_grid(lb = -Inf, rb = Inf, modes = modes_multi, f = f_multi, h = h_multi, h_prime = h_prime_multi)

plot(multi_grid)


multi_grid

multi_sampler = stors(multi_grid)

hist(multi_sampler(1000000))

save_grid(grid = multi_grid, grid_name = "multi_grid")

multi_grid = load_grid("multi_grid")

multi_sampler = stors(multi_grid)


hist(multi_sampler(1000000))

multi_trunc = trunc(multi_grid, -4,5)

hist(multi_trunc(1000000))


##  prebuilt dist

n = 10000


##    srnorm

srnorm_grid = grid_optimizer("srnorm" , verbose = TRUE)

print(srnorm_grid)

plot(srnorm_grid, x_min = -4 ,x_max =4 )

n = 10000

hist(srnorm(n))

ts = srnorm_truncate(-1,3)

trunc_sample = ts(n)

hist(trunc_sample)

min(trunc_sample)

max(trunc_sample)

library(RcppZiggurat)

m=100

microbenchmark::microbenchmark(
  rnorm(m),
  srnorm(m),
  zrnormR(m),
  times = 100
)

##    laplace

laplace_grid = grid_optimizer("srlaplace")

plot(laplace_grid)

hist(srlaplace(n))

tl = srlaplace_truncate(-2,3.3)

trunc_sample = tl(n)
  
hist(trunc_sample)


library(LaplacesDemon)

m=10000

microbenchmark::microbenchmark(
  srlaplace(m),
  LaplacesDemon::ralaplace(m),
  rexp(m),
  times = 100
)



##    srexp

srexp_grid = grid_optimizer("srexp")

plot(srexp_grid)

hist(srexp(n))

te=srexp_truncate(2,5)

trunc_sample = te(n)

hist(trunc_sample)

min(trunc_sample)

max(trunc_sample)

m=10000

microbenchmark::microbenchmark(
  srexp(m),
  rexp(m),
  times = 100
)





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
                       verbose = TRUE, target_sample_size = 100)

plot(norm_grid)

norm_grid

steps = 2000

norm_grid = build_grid(lb = -Inf, rb = Inf, modes , f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = steps)

norm_grid

plot(norm_grid)

norm_grid

norm_sampler = stors(norm_grid)

hist(norm_sampler(1000000))

norm_trunc = trunc(norm_grid, -2,2)

hist(norm_trunc(10000))

#=======================

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

plot(multi_grid, x_min = 0, x_max = 0.01)


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

srnorm_grid = grid_optimizer("srnorm", verbose = TRUE)

srnorm_grid = grid_optimizer("srnorm")

print(srnorm_grid)

plot(srnorm_grid, x_min = 3 ,x_max =7 )



plot(srnorm_grid,  x_min = 0 ,x_max= 0.01)

hist(srnorm(n))

trunc_srnorm = truncsrnorm(-1,3)

trunc_sample = trunc_srnorm(n)

hist(trunc_sample)

min(trunc_sample)

max(trunc_sample)

library(RcppZiggurat)

m=100000000

microbenchmark::microbenchmark(
  rnorm(m),
  srnorm(m),
  zrnormR(m),
  times = 1
)

##    srnorm

laplace_grid = grid_optimizer("laplace")

plot(laplace_grid)

plot(laplace_grid)

plot.grid

hist(laplace(n))

trunc_laplace = trunclaplace(-2,3.3)

trunc_sample = trunc_laplace(n)
  
hist(trunc_sample)


library(LaplacesDemon)

m=10000

microbenchmark::microbenchmark(
  laplace(m),
  LaplacesDemon::ralaplace(m),
  rexp(m),
  times = 100
)



##    srexp

srexp_grid = grid_optimizer("srexp")

plot(srexp_grid)

hist(srexp(n))

trunc_srexp=truncsrexp(2,5)

trunc_sample = trunc_srexp(n)

hist(trunc_sample)

min(trunc_sample)

max(trunc_sample)

m=10000

microbenchmark::microbenchmark(
  srexp(m),
  rexp(m),
  times = 100
)


##    srcauchy

sampe_size = 1

times = ceiling( 100000 / sampe_size)




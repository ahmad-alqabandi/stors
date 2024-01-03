
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

multi_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_multi, f = f_multi, h = h_multi, h_prime = h_prime_multi)

multi_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_multi, f = f_multi, a = 0.01, th = 0.001*2*2*2)
# multi_grid$grid_data$s_upper
# multi_grid$grid_data$p_a
plot(multi_grid)
multi_grid
# save_grid(multi_grid, "bi-modal")

# this is a grid with 000 infinity supported ... 
# srnorm(n, m, sd) call srsnorm 
# unit testing for grid(not sampling), micro in C (stencil in cales package),
# check for error when x from m_1 is greater than m_2


multi_sampler = stors(multi_grid)

hist(multi_sampler(1000000))

save_grid(grid = multi_grid, grid_name = "bimodal_3")

multi_grid = load_grid("multi_grid")

multi_sampler = stors(multi_grid)

hist(multi_sampler(1000000))



grid = stors::grid_optimizer("srnorm")

plot(grid)
# save_grid, load_grid

hist(srnorm(10000))


library(RcppZiggurat)

n=10000

microbenchmark::microbenchmark(
  rnorm(n),
  srnorm = srnorm(n),
  zrnormR(n),
  times = 100
)

-2.501206
-2.489405


microbenchmark::microbenchmark(
  zrnormR(1000),
  zrnormR(10000),
  zrnormR(100000),
  zrnormR(1000000),
  times = 100
)




vec =  grid$grid_data[1:(n),]$x
tail(vec, n = 3)
head(vec, n = 3)
length(vec)

vec = grid$grid_data[1:(n),]$s_upper
tail(vec, n = 3)
head(vec, n = 3)

length(vec)


vec =  grid$grid_data[2:(n-1),]$x
tail(vec, n = 3)
head(vec, n = 3)

length(vec)


vec =  grid$grid_data[5:(n-1),]$s_upper * grid$grid_data[5:(n-1),]$p_a
tail(vec, n = 3)
head(vec, n = 3)

length(vec)




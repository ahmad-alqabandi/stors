detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

detach_package(stors)

modes_multi = c(0,4)


p <- 0.5

q <- 0.5

f_multi <- function(x) {
  p * (sqrt(2 * pi))^(-1) * exp(-(x^2)/2) + q * ( sqrt(2 * pi))^(-1) * exp(-((x - 4)^2)/2)
}

h_multi <- function(x) {
  log(  p * (sqrt(2 * pi))^(-1) * exp(-(x^2)/2) + q * (sqrt(2 * pi))^(-1) * exp(-((x - 4)^2)/2))
}


h_prime_multi <- function(x) {
  (-(exp(-1/2 * (-4 + x)^2) * q * (-4 + x))/sqrt(2 * pi) - (exp(-x^2/2) * p * x)/sqrt(2 * pi))/((exp(-x^2/2) * p)/sqrt(2 * pi) + (exp(-1/2 * (-4 + x)^2) * q)/sqrt(2 * pi))
}

multi_grid = set_grid(lb = -Inf, rb = Inf, mode = modes_multi, f = f_multi, h = h_multi, h_prime = h_prime_multi)

multi_sampler = stors(multi_grid)

devtools::load_all()

# save_grid, load_grid

stors::grid_obtimizer("snorm")

hist(snorm(10000))


.Call(C_print_cached_grids)

multi_sampler = stors(multi_grid)


n=1
microbenchmark::microbenchmark(
  snorm(n),
  zrnormR(n),
  rnorm(n),
  times = 1000
)

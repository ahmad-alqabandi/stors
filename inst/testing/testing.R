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
  0.5 * (sqrt(2 * pi))^(-1) * exp(-(x^2)/2) + 0.5 * ( sqrt(2 * pi))^(-1) * exp(-((x - 4)^2)/2)
}

h_multi <- function(x) {
  log(  p * (sqrt(2 * pi))^(-1) * exp(-(x^2)/2) + q * (sqrt(2 * pi))^(-1) * exp(-((x - 4)^2)/2))
}


h_prime_multi <- function(x) {
  (-(exp(-1/2 * (-4 + x)^2) * q * (-4 + x))/sqrt(2 * pi) - (exp(-x^2/2) * p * x)/sqrt(2 * pi))/((exp(-x^2/2) * p)/sqrt(2 * pi) + (exp(-1/2 * (-4 + x)^2) * q)/sqrt(2 * pi))
}

multi_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_multi, f = f_multi, h = h_multi, h_prime = h_prime_multi)

multi_sampler = stors(multi_grid)

hist(multi_sampler(1000000))

save_grid(grid = multi_grid, grid_name = "bimodal_3")

multi_grid = load_grid("multi_grid")

multi_sampler = stors(multi_grid)

hist(multi_sampler(1000000))


stors::grid_obtimizer("srnorm")

hist(srnorm(n))

devtools::load_all()

# save_grid, load_grid


hist(srnorm(100000))


.Call(C_print_cached_grids)

multi_sampler = stors(multi_grid)


n=1000
microbenchmark::microbenchmark(
  srnorm(n),
  multi_sampler(n),
  zrnormR(n),
  rnorm(n),
  times = 100
)

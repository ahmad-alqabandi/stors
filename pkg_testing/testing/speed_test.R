
library(RcppZiggurat)
library(microbenchmark)

opt_cache_sizes = c(4 ,8 ,16 ,32, 64, 128, 256, 512, 1024)

opt_df_var = 4

opt_list_var = 20

opt_steps = round( ((opt_cache_sizes * 1024) - opt_list_var) / opt_df_var )

os = data.frame(steps = opt_steps, size = opt_cache_sizes)


m=10
times = 10000


# on.load opt steps = 65531, size = 256

microbenchmark::microbenchmark(
  srnorm(m),
  zrnormR(m),
  times = times
)


norm_grid = grid_optimizer("srnorm", steps = os$steps[2])


microbenchmark::microbenchmark(
  srnorm(m),
  zrnormR(m),
  times = times
)


norm_grid = grid_optimizer("srnorm", steps = os$steps[3])


microbenchmark::microbenchmark(
  srnorm(m),
  zrnormR(m),
  times = times
)


norm_grid = grid_optimizer("srnorm", steps = os$steps[4])


microbenchmark::microbenchmark(
  srnorm(m),
  zrnormR(m),
  times = times
)

norm_grid = grid_optimizer("srnorm", steps = os$steps[5])


microbenchmark::microbenchmark(
  srnorm(m),
  zrnormR(m),
  times = times
)

norm_grid = grid_optimizer("srnorm", steps = os$steps[6])


microbenchmark::microbenchmark(
  srnorm(m),
  zrnormR(m),
  times = times
)


norm_grid = grid_optimizer("srnorm", steps = os$steps[7])


microbenchmark::microbenchmark(
  srnorm(m),
  zrnormR(m),
  times = times
)


m <- 10^5

microbenchmark::microbenchmark(
  srbeta_custom(m),
  rbeta(m,shape1 = 2,shape2 = 3),
  times = 10
)


library(distributionsrd)

# Set parameters
n <- 1000  # Number of samples
shape <- 2      # Shape parameter
scale <- 1   # Scale parameter

# Generate random samples
samples <- rpareto(n, k = shape, xmin = scale)

plot(density(samples), main = "Density of Pareto Distributed Samples using distributionsrd Package")

grid <- srpareto_optimize(scale = scale, shape = shape, steps = 8000)


samples <- srpareto_custom(n)

plot(density(samples), main = "Density of Pareto Distributed Sample susing STORS Package")


microbenchmark(
  stors = srpareto_custom(n),
  rpareto = rpareto(n, k = shape, xmin = scale),
  times = 100
)
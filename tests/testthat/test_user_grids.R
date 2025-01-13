

test_that("MAIN", {

  modes_norm = 0
  f_norm <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x^2) }
  h_norm <- function(x) { log(f_norm(x)) }
  h_prime_norm <- function(x) { -x }

  norm_grid = build_grid(lb = -Inf, rb = Inf, mode = modes_norm,
                         f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)

  l_trunc <- 0
  u_trunc <- Inf
  norm_grid_trunc = build_grid(lb = l_trunc, rb = r_trunc, mode = modes_norm,
                         f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)


test_that("user smpling functions returns a sample of correct size and seeds", {


  sampling_function <- stors(norm_grid)

  set.seed(123)
  x1 <- sampling_function(10)

  set.seed(123)
  x2 <- sampling_function(10)


    expect_equal(x1, x2, info = "set.seed() does not works currectlly")


})


test_that("Truncated user grid samples test", {


  sampling_function_trunc <- stors(norm_grid_trunc)

  x <- sampling_function_trunc(1000)

  expect_true((min(x) >= l_trunc) && (max(x) <= u_trunc), info = paste0("generated samples must be within Normal distrebution\'s TRUNCATION BOUNDS"))

})



})

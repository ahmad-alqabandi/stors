

test_that("MAIN", {

  modes_norm <- 0
  f_norm <- function(x) 1 / sqrt(2 * pi) * exp(-0.5 * x^2)
  h_norm <- function(x) log(f_norm(x))
  h_prime_norm <- function(x) -x

  norm_grid <- build_proposal(lower = -Inf, upper = Inf, mode = modes_norm,
                              f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)

  l_trunc <- 0
  r_trunc <- Inf
  norm_grid_trunc <- build_proposal(lower = l_trunc, upper = r_trunc, mode = modes_norm,
                                    f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = 1000)


  test_that("user smpling functions returns a sample of correct size and seeds", {


    sampling_function <- build_sampler(norm_grid)

    set.seed(123)
    x1 <- sampling_function(10)

    set.seed(123)
    x2 <- sampling_function(10)


    expect_equal(x1, x2, info = "set.seed() does not works currectlly")


  })


  test_that("Truncated user grid samples test", {


    sampling_function_trunc <- build_sampler(norm_grid_trunc)

    x <- sampling_function_trunc(1000)

    expect_true((min(x) >= l_trunc) && (max(x) <= r_trunc), info = paste0("generated samples must be within Normal distrebution\'s TRUNCATION BOUNDS"))

  })



})

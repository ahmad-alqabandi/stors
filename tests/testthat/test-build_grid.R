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

# test_gris_d = grid_builder(a = 0.01, th = 0.4 ,mode = modes,f = f_norm, h = h_norm, h_prime = h_prime_norm )
# save(test_gris_d, file = "./tests/testthat/test_gris_d.Rdata")

# path = fs::path_package( "tests", "testthat", "test_gris_d.Rdata",package = "stors")

path <- testthat::test_path("test_gris_d.Rdata")
load(path)

test_that("multiplication works", {
  expect_equal(grid_builder(a = 0.01, th = 0.4, mode = modes, f = f_norm, h = h_norm, h_prime = h_prime_norm), test_gris_d)
})

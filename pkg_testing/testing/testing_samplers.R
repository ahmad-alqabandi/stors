
##    STORS uni modal


modes = 0

f_norm = function(x) {
   exp(-0.5 * x * x)
}

h_norm = function(x) {
  log(0.3989423) - (x * x) * (1 / 2)
}

h_prime_norm = function(x) {
  -x
}

norm_proposal <- build_proposal(lower = -Inf, upper = Inf, modes , f = f_norm,
                       steps = 1000, verbose = TRUE,  target_sample_size = 10000)

plot(norm_proposal)

s <- build_sampler(norm_proposal)

hist(s(10000))

? save_proposal

? print_proposals

? delete_proposal

? load_proposal

norm_proposal

save_proposal(norm_proposal, "normal")

print_proposals()

steps = 2000

norm_proposal = build_proposal(lower = -1, upper = Inf, modes , f = f_norm, h = h_norm, h_prime = h_prime_norm, steps = steps)

norm_proposal

plot(norm_proposal)

norm_proposal

norm_sampler = build_sampler(norm_proposal)

hist(norm_sampler(1000000))

norm_trunc = build_sampler(norm_proposal, -2,2)

hist(norm_trunc(10000))


##    STORS multi modal

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

multi_proposal = build_proposal(lower = -Inf, upper = Inf,verbose = T, mode = modes_multi, f = f_multi, h = h_multi, h_prime = h_prime_multi)

plot(multi_proposal)

multi_proposal

multi_sampler = build_sampler(multi_proposal)

hist(multi_sampler(1000000))

save_proposal(proposal = multi_proposal, proposal_name = "multi_proposal")

multi_proposal = load_proposal("multi_proposal")

multi_sampler = build_sampler(multi_proposal)

hist(multi_sampler(1000000))

delete_proposal("multi_proposal")


#======================= Built-in densities



##    srnorm

srnorm_proposal = proposal_optimizer("srnorm" , verbose = TRUE)

print(srnorm_proposal)

plot(srnorm_proposal, x_min = -4 ,x_max =4 )

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

laplace_proposal = proposal_optimizer("srlaplace")

plot(laplace_proposal)

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

microbenchmark::microbenchmark(
  srnorm_optimize(),
  times = 10
)

##    srexp

srexp_proposal = proposal_optimizer("srexp")

plot(srexp_proposal)

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



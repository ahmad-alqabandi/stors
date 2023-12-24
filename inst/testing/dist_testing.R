# Laplace distribution

b = 1
mu = 0
modes_Laplace = mu

f_Laplace <- function(x){
  
  1/(2*b) * exp(-abs(x - mu)/b)
  
}

h_Laplace <- function(x){log(f_Laplace(x))}

plot(seq(-100,100,0.01),h_Laplace(seq(-100,100,0.01)))

Laplace_grid = build_grid(mode = modes_Laplace , f = f_Laplace)

plot(Laplace_grid)

Laplace_sampler = stors(Laplace_grid)

hist(Laplace_sampler(100000) * 4, breaks = 100)

grid_optimizer("slaplace")


n = 1

microbenchmark::microbenchmark(
  Laplace_Stors = srlaplace(n),
  Laplace_C = stors::rlaplace_c(n),
  Laplace_ExtDist = ExtDist::rLaplace(n),
  Laplace_VGAM = VGAM::rlaplace(n),
  times = 1
)


hist(srlaplace(10000), breaks = 100)
hist(stors::rlaplace_c(10000), breaks = 100)


# Gumbel distribution

beta = 1
mu = 0
modes_gumbel = mu

f_gumbel <- function(x){
  
  exp(-(x + exp(-x)))
  
}

gumbel_grid = build_grid(mode = modes_gumbel , f = f_gumbel)

plot(gumbel_grid)

# Cauchy distribution
X0 = 0
gamma = 1
cauchy_mode = X0

f_cauchy <- function(x){
  
  1 / (pi*(1+x^2))
  
}

plot(seq(-100,100,0.01),f_cauchy(seq(-100,100,0.01)), 'l')

h_cauchy <- function(x) {log(f_cauchy(x))}

plot(seq(-100,100,0.01),h_cauchy(seq(-100,100,0.01)), 'l')

h_prime_cauchy <- function(x){
  
  -(2*x)/(1+x^2)
  
}

cauchy_grid = build_grid(mode = cauchy_mode , f = f_cauchy, th = 0.01)

plot(cauchy_grid)

cauchy_grid_hp = build_grid(mode = cauchy_mode , f = f_cauchy,h = h_cauchy, h_prime = h_prime_cauchy, th = 0.01)

plot(cauchy_grid_hp)



# test truncation

xl = -5
huxl = normg$lt_properties[5] * ( xl -  normg$grid_data$x[1]) + normg$lt_properties[3]

uxl =( (normg$lt_properties[4]) * (exp(huxl) - normg$lt_properties[1]) ) / (sum(normg$areas))

normg = grid_optimizer("srnorm")

hist(srnorm(100000))

hist(truncsrnorm(-2,2)(100000))

trnc23=truncsrnorm(-2,3)
trnc23_paste0 =truncsrnorm_paste_0(-2,3)
trnc23_hard_code = truncsrnorm_hard_code(-2,3)
trnc23_trick =  truncsrnorm_trick(-2,3)


n = 1
library(truncnorm)
microbenchmark::microbenchmark(
  trnc23_trick(n),
  rtruncnorm(n, a=-2, b=3, mean = 0, sd = 1),
  trnc23(n),
  trnc23_paste0(n),
  trnc23_hard_code(n),
  truncsrnorm_inter(n,0.02325644 ,0.99846791),
  times = 1000
)

microbenchmark::microbenchmark(
  truncsrnorm_trick(-2,3)(n),
  rtruncnorm(n, a=-2, b=3, mean = 0, sd = 1),
  truncsrnorm(-2,3)(n),
  truncsrnorm_paste_0(-2,3)(n),
  truncsrnorm_hard_code(-2,3)(n),
  truncsrnorm_inter(n,0.02325644 ,0.99846791),
  times = 1000
)

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



#
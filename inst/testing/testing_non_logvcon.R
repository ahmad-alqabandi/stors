

##    laplace


f = function(x) 0.5 * exp(-abs(x))

h = function(x) -0.6931472 -abs(x)

h_prime = function(x) {
  if(x > 0){
    -1
  }else{
    1
  }
} 

laplace_grid = build_grid(lb = -Inf, rb = Inf, modes = 0 , f = f, h = h, h_prime = h_prime,
                       steps = 1000, verbose = TRUE,  target_sample_size = 10000)

plot(laplace_grid)

laplace_sampler = stors(laplace_grid)

hist(laplace_sampler(100000))

plot(laplace_grid, -5.5, -5.49)

x1 = laplace_grid$grid_data$x[1]

h_prime(x1)

g =  function (x) exp( h_prime(x1) * ( x - x1 ) + h(x1))


l_range = c(-Inf, x1)

g_l =  function (x) exp( x + -0.693147 )

f_l =  function (x) 0.5 * exp(x)

w = function(x) g_l(x) - f_l(x)

w = function(x) exp( x + -0.693147 ) -  0.5 * exp(x)

w(l_range[1])
w(l_range[2])

xx = seq(from = -10, to = l_range[2], length.out = 1000)

plot(xx, w(xx), type = "line")

h_upper <- function(grid_point, val, h_prime, h) {
  h_prime(grid_point) * (val - grid_point) + h(grid_point)
}




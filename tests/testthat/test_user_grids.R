


dists <- list(
  unimolad_normal = list(
    lb = -Inf,
    rb = Inf,
    modes = 0,
    f = function(x) {
      0.3989423 * exp(-0.5 * x * x)
    }
  ),
  multimodal_normal = list(
    lb = -Inf,
    rb = Inf,
    modes = c(-2, 2),
    f = function(x) {
      w1 <- 0.5
      w2 <- 0.5
      sd1 <- 0.5
      sd2 <- 0.5
      mean1 <- -2
      mean2 <- 2
      w1 * (1 / (sqrt(2 * pi) * sd1)) * exp(-((x - mean1)^2) / (2 * sd1^2)) +
        w2 * (1 / (sqrt(2 * pi) * sd2)) * exp(-((x - mean2)^2) / (2 * sd2^2))
    }
  )
)


steps = 2048

grids <- list()

for( name in names(dists) ){
  
  grids[[name]] <- do.call(build_grid, c(dists[[name]],list(steps = steps)))
  
}


test_that("User_grid smpling functions returns a sample of correct size and seeds",{
  
  
  for( name in names(grids)){

    sampler_fun <- stors(grids[[name]])
    x <- sampler_fun(10)
    n <- length(x)
    
    expect_equal(n,10, info =paste0(name," returned sample size does not match sample size arg"))
    
    set.seed(1234)
    x1 <- sampler_fun(10)
    set.seed(1234)
    x2 <- sampler_fun(10)
    
    expect_equal(x1,x2, info = "set.seed() does not works as intendded")
    
  }  
  
})


test_that("User_grid sampling functions, samples properties tests", {
  for (name in names(grids)) {
    
    lb <- grids[[name]]$grid_bounds[1]
    rb <- grids[[name]]$grid_bounds[2]
    
    sampler_fun <- stors(grids[[name]])
    x <- sampler_fun(1000)
    
    expect_true((min(x) >= lb) && (max(x) <= rb) , info = 'generated samples must be within the distrebution\'s BOUNDS')
    
    
    rnds <- runif(2)
    rnds <- rnds[order(rnds)]
    
    poss_min <- if(is.infinite(lb)) -10 else lb
    poss_max <- if(is.infinite(rb)) 10 else rb
    
    l_trunc <- poss_min + (poss_max - poss_min) * rnds[1]
    u_trunc <- poss_min + (poss_max - poss_min) * rnds[2]
    

    trunc_sampler_fun <- stors(grids[[name]], l_trunc, u_trunc)
    x <- trunc_sampler_fun(10000)
    
    expect_true((min(x) >= l_trunc) && (max(x) <= u_trunc) ,
                info = paste0('generated samples must be within ', name,
                              ' distrebution\'s TRUNCATION BOUNDS'))
  }
  
})
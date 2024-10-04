# BG built-in grids

test_that("Built-in smpling functions returns a sample of correct size",{
  
  for( name in stors:::stors_env$grids$builtin$names){
    
    srname_txt <- paste0(name,'(10)')
    srname_exp <- parse(text = srname_txt)
    x <- eval(srname_exp)
    n <- length(x)
    
    expect_equal(n,10, info =paste0(name," returned sample size does not match sample size arg"))
    
  }
  
  
})

test_that("BG: srnorm responed correcttly to `set.seed()` ",{
  
  set.seed(1234)
  x1 <- srnorm(3)
  set.seed(1234)
  x2 <- srnorm(3)
  
  expect_equal(x1,x2)
  
})



test_that("Built-in sampling functions, samples properties tests", {
  for (name in stors:::stors_env$grids$builtin$names) {

    lb <- pbgrids[[name]]$lb
    rb <- pbgrids[[name]]$rb
    
    srname_txt <- paste0(name, '(1000)')
    srname_exp <- parse(text = srname_txt)
    x <- eval(srname_exp)
    
    expect_true((min(x) >= lb) && (max(x) <= rb) , info = 'generated samples must be within the distrebution\'s BOUNDS')
    
    
    rnds <- runif(2)
    rnds <- rnds[order(rnds)]
    
    poss_min <- ifelse(is.finite(lb), lb, pbgrids[[name]]$is.lb)
    poss_max <- ifelse(is.finite(rb), rb, pbgrids[[name]]$is.rb)
    
    l_trunc <- poss_min + (poss_max - poss_min) * rnds[1]
    u_trunc <- poss_min + (poss_max - poss_min) * rnds[2]
    
    srname_truncated_fun_txt <- paste0(name, '_truncate(',l_trunc,',',u_trunc,')','(10000)')
    srname_truncated_fun_exp <- parse(text = srname_truncated_fun_txt)
    x <- eval(srname_truncated_fun_exp)
    
    expect_true((min(x) >= l_trunc) && (max(x) <= u_trunc) ,
                info = paste0('generated samples must be within ', name,
                              ' distrebution\'s TRUNCATION BOUNDS'))
    
  }
  
})





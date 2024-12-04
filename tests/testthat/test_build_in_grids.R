# BG built-in grids

test_that("Built_in smpling functions returns a sample of correct size and seeds",{
  
  for( name in names(stors:::pbgrids)){
    srname_txt <- paste0(name,'(10)')
    srname_exp <- parse(text = srname_txt)
    x <- eval(srname_exp)
    n <- length(x)
    
    expect_equal(n,10, info =paste0(name," returned sample size does not match sample size arg"))
    
    set.seed(1234)
    x1 <- eval(srname_exp)
    set.seed(1234)
    x2 <- eval(srname_exp)
    
    expect_equal(x1,x2, info = "set.seed() does not works currectlly")
  }
  
})


test_that("Built_in sampling functions, samples properties tests", {
  for (name in names(stors:::pbgrids)) {

    lb <- stors:::pbgrids[[name]]$lb
    rb <- stors:::pbgrids[[name]]$rb
    
    rnds <- runif(2)
    rnds <- rnds[order(rnds)]
    
    poss_min <- if(is.infinite(lb)) -10 else lb
    poss_max <- if(is.infinite(rb)) 10 else rb
    
    l_trunc <- poss_min + (poss_max - poss_min) * rnds[1]
    u_trunc <- poss_min + (poss_max - poss_min) * rnds[2]
    
    
    g <- eval(parse(text = paste0(name,"_optimize( xl = ",l_trunc,",xr = ",u_trunc,")")))
    
    
    x <- eval(parse(text = paste0(name, '(1000)')))
    
    expect_true((min(x) >= l_trunc) && (max(x) <= u_trunc) , info = paste0('generated samples must be within ',
                                                                           name, ' distrebution\'s TRUNCATION BOUNDS'))    
    # 
    # rnds <- runif(2)
    # rnds <- rnds[order(rnds)]
    # 
    # poss_min <- if(is.infinite(lb)) -10 else lb
    # poss_max <- if(is.infinite(rb)) 10 else rb
    # 
    # l_trunc <- poss_min + (poss_max - poss_min) * rnds[1]
    # u_trunc <- poss_min + (poss_max - poss_min) * rnds[2]
    # 
    # srname_truncated_fun_txt <- paste0(name, '_truncate(',l_trunc,',',u_trunc,')','(10000)')
    # srname_truncated_fun_exp <- parse(text = srname_truncated_fun_txt)
    # x <- eval(srname_truncated_fun_exp)
    # 
    # expect_true((min(x) >= l_trunc) && (max(x) <= u_trunc) ,
    #             info = paste0('generated samples must be within ', name,
    #                           ' distrebution\'s TRUNCATION BOUNDS'))
  }
  
})

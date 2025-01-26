# BG built-in grids

test_that("Built_in smpling functions returns a sample of correct size and seeds", {

  for (name in names(built_in_proposals)) {

    if (built_in_proposals[[name]]$scalable) {
      srname_txt <- paste0(name, "(10)")
    } else {
      srname_txt <- paste0(name, "_custom(10)")
    }

    srname_exp <- parse(text = srname_txt)
    x <- eval(srname_exp)
    n <- length(x)

    expect_equal(n, 10, info = paste0(name, " returned sample size does not match sample size arg"))

    set.seed(1234)
    x1 <- eval(srname_exp)
    set.seed(1234)
    x2 <- eval(srname_exp)

    expect_equal(x1, x2, info = "set.seed() does not works currectlly")
  }

})


test_that("Built_in sampling functions, samples properties tests", {
  for (name in setdiff(names(built_in_proposals), "srpareto")) {

    lb <- built_in_proposals[[name]]$lb
    rb <- built_in_proposals[[name]]$rb

    rnds <- runif(2)
    rnds <- rnds[order(rnds)]

    poss_min <- if (is.infinite(lb)) -10 else lb
    poss_max <- if (is.infinite(rb)) 10 else rb

    l_trunc <- poss_min + (poss_max - poss_min) * rnds[1]
    u_trunc <- poss_min + (poss_max - poss_min) * rnds[2]


    g <- eval(parse(text = paste0(name, "_optimize( xl = ", l_trunc, ",xr = ", u_trunc, ")")))


    if (built_in_proposals[[name]]$scalable) {
      srname_txt <- paste0(name, "(1000)")
    } else {
      srname_txt <- paste0(name, "_custom(1000)")
    }

    srname_exp <- parse(text = srname_txt)
    x <- eval(srname_exp)

    expect_true((min(x) >= l_trunc) && (max(x) <= u_trunc), info = paste0("generated samples must be within ",
                                                                           name, " distrebution\'s TRUNCATION BOUNDS"))
  }

  for (name in names(built_in_proposals)) {

    fun_name <- paste0(name, "_optimize")
    do.call(fun_name, list(steps = 4091))

  }

})

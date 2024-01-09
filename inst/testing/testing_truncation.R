
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

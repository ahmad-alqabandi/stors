sample_size <- 1000

shape <- 1
scale <- .5

srgamma_optimize(shape = shape, scale = scale)

df <- data.frame( x = c(rgamma(sample_size, shape = shape, scale = scale),
                          srgamma(sample_size)) ,
                   method = c(rep("BASE",sample_size),rep("STORS",sample_size)))


ggplot(data = df, aes(x = x, colour = method)) +
  geom_histogram(alpha = 0.2, position = "identity", binwidth = 0.5)


microbenchmark::microbenchmark(
  srgamma(sample_size),
  rgamma(sample_size, shape = shape, scale = scale),
  times = 1000
)

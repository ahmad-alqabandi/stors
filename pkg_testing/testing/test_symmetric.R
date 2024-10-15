# library(ggplot2)
# 
# n <- 11
# 
# steps <- 2^(1:n)
# 
# sample_size = 100000
# 
# df <- data.frame(srnorm = rep(NA, n), srnorm_symmetric = rep(NA, n), steps = rep(NA, n))
# 
# for( i in (1:length(steps))){
# 
# norm <- srnorm_optimize(steps = steps[i])
# 
# norm_symmetric <- srnorm_symmetric_optimize( steps = steps[i] * 2)
# 
# res_df <- microbenchmark::microbenchmark(
#   norm = srnorm(sample_size),
#   norm_symmetric = srnorm_symmetric(sample_size),
#   times = 100
# )
# 
# mu_norm <- mean(res_df[res_df$expr == "norm", "time"])
# 
# mu_symmetric <- mean(res_df[res_df$expr == "norm_symmetric", "time"])
# 
# df[i,] <- c(mu_norm,mu_symmetric, steps[i])
# 
# }
# 
# df_long <- data.frame(
#   steps = rep(df$steps, 2),
#   value = c(df$srnorm, df$srnorm_symmetric),
#   type = rep(c("srnorm", "srnorm_symmetric"), each = length(df$steps))
# )
# 
# p1 <- ggplot(data = df_long, aes(x = steps, y = value, color = type)) +
#   geom_line() +
#   scale_color_manual(values = c("srnorm" = "red", "srnorm_symmetric" = "green")) +
#   # scale_x_continuous(
#   #   breaks = df$steps,
#   #   labels = df$steps
#   # ) +
#   scale_x_log10(breaks = df$steps) +
#   scale_y_log10() +
#   labs(x = "number of steps", y = "average time", color = "sampling method")
# 
# print(p1)




# ==============================================================================



library(RcppZiggurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(microbenchmark)
library(dplyr)
library(cowplot)

sample_size <- 10000

srnorm_optimize()

old_srnorm_optimize()


microbenchmark::microbenchmark(
  srnorm(sample_size),
  old_srnorm(sample_size),
  zrnormR(sample_size),
  times = 1000
)

mu =2; sd = 1.5

srnorm_optimize(mu, sd, steps = 4091)

sample_size %>% srnorm() %>% hist()
sample_size %>% old_srnorm_scaled( mu, sd) %>% hist()

microbenchmark::microbenchmark(
  srnorm(sample_size),
  old_srnorm_scaled(sample_size, mu, sd),
  times = 1000
)



df <- 5

srchisq_optimize(df)

data1 <- data.frame(value = rchisq(sample_size, df))
data2 <- data.frame(value = srchisq(sample_size))

p1 <- ggplot(data1, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black", alpha = 0.7) +
  ggtitle(paste0("rchisq, df = ",df)) +
  theme_minimal()

p2 <- ggplot(data2, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "salmon", color = "black", alpha = 0.7) +
  ggtitle(paste0("STORS srchisq, df = ",df)) +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)


res <- microbenchmark::microbenchmark(
  srchisq(sample_size),
  rchisq(sample_size,df),
  times = 1000
)

levels(res$expr) <- c("STORS srchisq", "STATS rchisq")

# Convert benchmark results to a data frame
benchmark_df <- as.data.frame(res)

# Summarize benchmark results
benchmark_summary <- benchmark_df %>%
  group_by(expr) %>%
  summarise(
    min_time = min(time) / 1e6,    # Convert nanoseconds to milliseconds
    median_time = median(time) / 1e6,
    max_time = max(time) / 1e6
  ) %>%
  rename(
    Method = expr,
    Min_ms = min_time,
    Median_ms = median_time,
    Max_ms = max_time
  )

# Round the time values for better display
benchmark_summary <- benchmark_summary %>%
  mutate(
    Min_ms = round(Min_ms, 3),
    Median_ms = round(Median_ms, 3),
    Max_ms = round(Max_ms, 3)
  )

# Create table grob
benchmark_table <- tableGrob(benchmark_summary)

# Arrange the histograms side by side
histograms <- arrangeGrob(p1, p2, ncol = 2)

# Arrange the histograms and the table
combined_plot <- grid.arrange(
  histograms,
  benchmark_table,
  nrow = 2,
  heights = c(2, 1)  # Adjust heights as needed
)

# Display the combined plot
print(combined_plot)




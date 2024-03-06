rm(list = ls())

library(ggplot2)
library(dplyr)

source("figures_scripts/illustration_grhmc_median_crossing_times.R")

# consider a N(0, Sigma) target where
Sigma <- matrix(
  c(10, 5,
    5, 10^3),
  nrow = 2,
  byrow = TRUE
)

mu <- c(0, 0)

prec <- solve(Sigma)

# target gradient (change to get another target distribution)
target.grad <- function(q) {return(-prec%*%(q - mu))}

# Use desired time to median as pi for now. 
# Also: Use lambda = 0.3
lambda <- 0.3
time_period <- 10000
n_samples <- 100000

original_sim_run <- illustration_grhmc_median_crossing_times(
  log_target_grad = target.grad,
  lambda = lambda,
  T = time_period,
  n_samples = n_samples,
  median_elements = c(0, 0),
  s_elements = c(sqrt(10), sqrt(1000)),
  random_state = 42
)

original_sim_output <- original_sim_run$output_from_ode_solver

quantile(diff(original_sim_run$time_q_root_median[[1]]))
quantile(diff(original_sim_run$time_q_root_median[[2]]))

summary(original_sim_run$q_samples)
cov(original_sim_run$q_samples)

plot(original_sim_output$time[original_sim_output$time <= 500], original_sim_run$q_samples[original_sim_output$time <= 500, 1])
abline(h = 0, col = "red")

plot(original_sim_output$time[original_sim_output$time <= 500], original_sim_run$q_samples[original_sim_output$time <= 500, 2])
abline(h = 0, col = "red")

# ggplot2 version

# check the first 500 time units - corresponds to Figure 1

df_first_coordinate <- data.frame(time = original_sim_output$time[original_sim_output$time <= 500], q1 = original_sim_run$q_samples[original_sim_output$time <= 500, 1])
df_second_coordinate <- data.frame(time = original_sim_output$time[original_sim_output$time <= 500], q2 = original_sim_run$q_samples[original_sim_output$time <= 500, 2])

df_first_coordinate %>%
  ggplot(aes(x = time, y = q1)) +
  geom_line() + 
  geom_hline(yintercept = 0, col = "red") + 
  xlab("Time") + 
  ylab(bquote(q[1])) +
  # ggtitle(expression(paste("First coordinate: ", mu, " = 0, ", sigma^2, " = 10", sep = ""))) +
  # theme(plot.title = element_text(hjust = 0.5))
  theme(
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 32)
  )

df_second_coordinate %>%
  ggplot(aes(x = time, y = q2)) +
  geom_line() + 
  geom_hline(yintercept = 0, col = "red") + 
  xlab("Time") + 
  ylab(bquote(q[2])) + 
  # ggtitle(expression(paste("Second coordinate: ", mu, " = 0, ", sigma^2, " = 100", sep = ""))) +
  # theme(plot.title = element_text(hjust = 0.5))
  theme(
    axis.text = element_text(size = 28),
    axis.title = element_text(size = 32)
  )


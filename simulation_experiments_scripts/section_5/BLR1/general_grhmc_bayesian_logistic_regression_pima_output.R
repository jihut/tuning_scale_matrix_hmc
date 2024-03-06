rm(list = ls())

library(MASS)
source("simulation_experiments_scripts/general_scripts/general_grhmc_rstan_monitor_transformed_samples_stored.R")
source("simulation_experiments_scripts/general_scripts/single_model_rstan_summary_table.R")

data("Pima.te")
data("Pima.tr")

Pima_full <- rbind(Pima.te, Pima.tr)

head(Pima_full)

Pima_full[, 1:7] <- scale(Pima_full[, 1:7]) # scale the numeric predictor variables

Pima_full[, "intercept"] <- rep(1, nrow(Pima_full))

x_matrix <- as.matrix(Pima_full[, -which(colnames(Pima_full) == "type")])
y_vector <- as.integer(Pima_full$type == "Yes")

transposed_x_matrix <- t(x_matrix) 

alpha <- 100 # sigma^2, between 10^2 and 100^2 usually

log_target_gradient_with_gaussian_prior <- function(beta){
  
  eta <- x_matrix %*% beta
  
  p <- 1 / (1 + exp(-eta))
  
  likelihood_term <- transposed_x_matrix %*% (y_vector - p)
  
  log_prior <- - beta / alpha
  
  return((likelihood_term + log_prior))
  
}

general_grhmc_bayesian_logistic_regression_pima_output <- general_grhmc_rstan_monitor_samples_stored(
  export_result = "simulation_experiments_scripts/section_5/BLR1/general_grhmc_bayesian_logistic_regression_pima_output.RDS",
  log_target_grad = log_target_gradient_with_gaussian_prior,
  lambda_initial = 0.2,
  n_chains = 10,
  n_generated_samples = 100000,
  time_period_generating_samples = 100000,
  n_parameters = 8,
  desired_time_until_crossing_median = pi,
  random_state = 1,
  vari_and_isg_time_period_adaptive = 6000,
  mct_step_sim_time_first_run = 1000,
  mct_step_sim_time_second_run = 5000,
  mct_step_gamma_elements_dual_averaging_first_run = rep(10, 8), 
  mct_step_gamma_elements_dual_averaging_second_run = rep(25, 8),
  mct_step_mu_log_scale_second_run = TRUE,
  mct_step_mu_scale_factor_second_run = 1.1,
  lambda_adaptive_step_sim_time = 5000
)

general_grhmc_bayesian_logistic_regression_pima_output <- readRDS("simulation_experiments_scripts/section_5/BLR1/general_grhmc_bayesian_logistic_regression_pima_output.RDS") # If the above has already been run, then read the generated object using this line

# RStan summary table where 50000 out of 100000 samples are considered from each chain

rstan_table <- rstan_summary_table(RDS_object = general_grhmc_bayesian_logistic_regression_pima_output, sample_grid = 2, name = "BLR1")
rstan_table

saveRDS(rstan_table$min_max_summary_table, file = "simulation_experiments_scripts/section_5/BLR1/relevant_rstan_table_BLR1.RDS")

# Check correlation between the parameters using samples from VARI

samples_VARI_df <- data.frame(init = numeric(500000)) 

for(i in 1:8){
  
  name <- paste("VARI_nr", i, sep = "")
  
  samples_VARI_df[name] <- as.numeric(general_grhmc_bayesian_logistic_regression_pima_output$VARI_3d_array[seq(from = 2, to = 100000, by = 2), , i])
  
  if(i == 8){
    samples_VARI_df <- samples_VARI_df[, -(colnames(samples_VARI_df) %in% c("init"))]
  }
  
}

head(samples_VARI_df)

sort(abs(unique(as.numeric(cor(samples_VARI_df)))), decreasing = T)

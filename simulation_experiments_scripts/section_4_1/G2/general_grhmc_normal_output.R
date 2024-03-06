rm(list = ls())

source("simulation_experiments_scripts/general_scripts/general_grhmc_rstan_monitor_transformed_samples_stored.R")
source("simulation_experiments_scripts/general_scripts/single_model_rstan_summary_table.R")

# Bivariate normal distribution

source("simulation_experiments_scripts/section_4_1/G2/normal_log_target_grad.R")

general_grhmc_normal_output <- general_grhmc_rstan_monitor_samples_stored(
  export_result = "simulation_experiments_scripts/section_4_1/G2/general_grhmc_normal_output.RDS",
  log_target_grad = target.grad,
  lambda_initial = 0.2,
  n_chains = 10,
  n_generated_samples = 100000,
  time_period_generating_samples = 100000,
  n_parameters = 2,
  desired_time_until_crossing_median = pi,
  random_state = 1,
  vari_and_isg_time_period_adaptive = 6000,
  mct_step_sim_time_first_run = 1000,
  mct_step_sim_time_second_run = 5000,
  mct_step_gamma_elements_dual_averaging_first_run = rep(10, 2), 
  mct_step_gamma_elements_dual_averaging_second_run = rep(25, 2),
  mct_step_mu_log_scale_second_run = TRUE,
  mct_step_mu_scale_factor_second_run = 1.1,
  lambda_adaptive_step_sim_time = 5000
)

general_grhmc_normal_output <- readRDS("simulation_experiments_scripts/section_4_1/G2/general_grhmc_normal_output.RDS") # If the above has already been run, then read the generated object using this line

# RStan summary table where 50000 out of 100000 samples are considered from each chain

rstan_table <- rstan_summary_table(RDS_object = general_grhmc_normal_output, sample_grid = 2, name = "G2")
rstan_table

saveRDS(rstan_table$all_summary_table, file = "simulation_experiments_scripts/section_4_1/G2/relevant_rstan_table_G2.RDS")

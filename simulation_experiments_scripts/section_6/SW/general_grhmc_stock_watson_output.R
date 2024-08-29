rm(list = ls())

source("simulation_experiments_scripts/general_scripts/ver2_general_grhmc_rstan_monitor_transformed_samples_stored.R")
source("simulation_experiments_scripts/general_scripts/single_model_rstan_summary_table.R")

# Stock and Watson model 

library(bridgestan)

y <- as.vector(read.table("simulation_experiments_scripts/section_6/SW/USdata_updated.txt")$x)

T <- length(y)

jsonlite::write_json(list(T = T, y = y), "simulation_experiments_scripts/section_6/SW/data.json")

#cc <- stanc_builder(file="sw_drhmc_tauc.stan",allow_undefined=TRUE, isystem=CIP_header_path())
#drhmc.mdlt <- stan_model(stanc_ret=cc,allow_undefined=TRUE,includes=CIP_include())

sw_model <- bridgestan::StanModel$new("simulation_experiments_scripts/section_6/SW/sw_drhmc_tauc.stan", "simulation_experiments_scripts/section_6/SW/data.json", 1234) # use BridgeStan to obtain the gradient of log posterior of the model
n_parameters <- 1 + T - 1 + T + T # based on stan file
sw_model$log_density_gradient(c(rep(0, 1 + T - 1 + T + T)))

target.grad <- function(q) {
  
  sw_model$log_density_gradient(q, propto = FALSE, jacobian = FALSE)$gradient
  
}

general_grhmc_stock_watson_output <- general_grhmc_rstan_monitor_samples_stored(
  export_result = "simulation_experiments_scripts/section_6/SW/general_grhmc_stock_watson_output.RDS",
  log_target_grad = target.grad,
  lambda_initial = 0.2,
  n_chains = 10,
  n_generated_samples = 20000,
  time_period_generating_samples = 10000,
  n_parameters = n_parameters,
  desired_time_until_crossing_median = pi,
  random_state = 1,
  vari_and_isg_time_period_adaptive = 6000,
  mct_step_sim_time_first_run = 1000,
  mct_step_sim_time_second_run = 5000,
  mct_step_gamma_elements_dual_averaging_first_run = rep(10, n_parameters), 
  mct_step_gamma_elements_dual_averaging_second_run = rep(25, n_parameters),
  mct_step_mu_log_scale_second_run = TRUE,
  mct_step_mu_scale_factor_second_run = 1.1,
  mct_step_time_step_update_data_median = 10,
  lambda_adaptive_step_sim_time = 1000
)

general_grhmc_stock_watson_output <- readRDS("simulation_experiments_scripts/section_6/SW/general_grhmc_stock_watson_output.RDS") # If the above has already been run, then read the generated object using this line

# RStan summary table where 10000 out of 20000 samples are considered from each chain

rstan_table <- rstan_summary_table(RDS_object = general_grhmc_stock_watson_output, sample_grid = 2, name = "SW")
rstan_table

saveRDS(rstan_table, file = "simulation_experiments_scripts/section_6/SW/relevant_rstan_table_SW.RDS")


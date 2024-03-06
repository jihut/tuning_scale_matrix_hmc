# Entire MCT procedure
# First step: Tuning m and S
# Second step: Tuning lambda
# Third step: Generating samples based on fixed m, S and lambda from the previous two steps

source("tuning_procedures_main_scripts/MCT/first_step_MCT_general_function_of_grhmc_adaptive.R")
source("tuning_procedures_main_scripts/general_scripts/second_step_lambda_adapt_function_lsodar_of_grhmc_adaptive.R")
source("tuning_procedures_main_scripts/general_scripts/third_step_general_function_of_grhmc_adaptive.R")

final_MCT_general_function_of_grhmc_adaptive <- function(
  log_target_grad,
  lambda_initial = 0.2,
  n_parameters = NULL,
  median_crossing_step_sim_time_first_run = 1000,
  median_crossing_step_sim_time_second_run = 5000,
  median_crossing_step_n_samples = 1000,
  desired_time_until_crossing_median = pi,
  median_crossing_step_diag_s_elements_initial = NULL,
  median_crossing_step_m_initial = NULL,
  median_crossing_step_qbar_initial = NULL,
  median_crossing_step_pbar_initial = NULL,
  random_state = NULL,
  maxsteps = NULL,
  rtol = NULL,
  atol = NULL,
  median_crossing_step_proportion_time_until_adaptive = 0.05,
  median_crossing_step_time_step_update_data_median = 1,
  median_crossing_step_kappa_elements_dual_averaging = NULL,
  median_crossing_step_gamma_elements_dual_averaging_first_run = NULL,
  median_crossing_step_gamma_elements_dual_averaging_second_run = NULL,
  median_crossing_step_t0_elements_dual_averaging = NULL,
  median_crossing_step_mu_elements_dual_averaging_first_run = NULL,
  median_crossing_step_mu_log_scale_second_run = TRUE,
  median_crossing_step_mu_scale_factor_second_run = NULL,
  lambda_adaptive_step_sim_time = 5000,
  lambda_adaptive_step_ema_beta = 0.99,
  lambda_adaptive_step_factor = 1,
  lambda_adaptive_step_n_samples = 1000,
  lambda_adaptive_step_max_nut_time = 100,
  final_sampling_step_sim_time = 5000,
  final_sampling_step_n_samples = 5000
  
){
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  first_step_run <- first_step_grhmc_transformed_function_MCT(
    log_target_grad = log_target_grad,
    lambda = lambda_initial,
    n_parameters = n_parameters,
    sim_time_first_run = median_crossing_step_sim_time_first_run,
    sim_time_second_run = median_crossing_step_sim_time_second_run,
    n_samples = median_crossing_step_n_samples,
    desired_time_until_crossing_median = desired_time_until_crossing_median,
    diag_s_elements_initial = median_crossing_step_diag_s_elements_initial,
    m_initial = median_crossing_step_m_initial,
    qbar_initial = median_crossing_step_qbar_initial,
    pbar_initial = median_crossing_step_pbar_initial,
    maxsteps = maxsteps,
    rtol = rtol,
    atol = atol,
    proportion_time_until_adaptive_start = median_crossing_step_proportion_time_until_adaptive,
    time_step_update_data_median = median_crossing_step_time_step_update_data_median,
    kappa_elements_dual_averaging = median_crossing_step_kappa_elements_dual_averaging,
    gamma_elements_dual_averaging_first_run = median_crossing_step_gamma_elements_dual_averaging_first_run,
    gamma_elements_dual_averaging_second_run = median_crossing_step_gamma_elements_dual_averaging_second_run,
    t0_elements_dual_averaging = median_crossing_step_t0_elements_dual_averaging,
    mu_elements_dual_averaging_first_run = median_crossing_step_mu_elements_dual_averaging_first_run,
    mu_log_scale_second_run = median_crossing_step_mu_log_scale_second_run,
    mu_scale_factor_second_run = median_crossing_step_mu_scale_factor_second_run
  )  
  
  second_step_run <- second_step_lambda_adapt_function(
    fun = log_target_grad,
    q0 = first_step_run$second_median_crossing_run$qbar_end,
    p0 = first_step_run$second_median_crossing_run$pbar_end,
    m.c = first_step_run$second_median_crossing_run$m_adaptive,
    S.c = first_step_run$second_median_crossing_run$s_adaptive,
    Tmax = lambda_adaptive_step_sim_time,
    lambda0 = lambda_initial,
    ema.beta = lambda_adaptive_step_ema_beta,
    fac = lambda_adaptive_step_factor,
    nsample = lambda_adaptive_step_n_samples,
    max.nut.time = lambda_adaptive_step_max_nut_time,
    maxsteps = maxsteps,
    rtol = rtol,
    atol = atol
  )
  
  final_step_run <- third_step_grhmc_transformed_function(
    log_target_grad = log_target_grad,
    lambda = second_step_run$lambda.last,
    T = final_sampling_step_sim_time,
    n_samples = final_sampling_step_n_samples,
    diag_s_elements_initial = first_step_run$second_median_crossing_run$s_adaptive,
    m_initial = first_step_run$second_median_crossing_run$m_adaptive,
    qbar_initial = second_step_run$q.last,
    pbar_initial = second_step_run$p.last,
    Lambda_initial = 0,
    u_initial = rexp(1),
    maxsteps = maxsteps,
    rtol = rtol,
    atol = atol
  )
    
  
  return(
    list(
      first_step_run = first_step_run,
      second_step_run = second_step_run,
      final_step_run = final_step_run
    )
  )

}

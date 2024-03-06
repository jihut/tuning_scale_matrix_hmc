# Entire VARI procedure
# First step: Tuning m and S
# Second step: Tuning lambda
# Third step: Generating samples based on fixed m, S and lambda from the previous two steps

source("tuning_procedures_main_scripts/VARI/first_step_VARI_general_function_of_grhmc_adaptive.R")
source("tuning_procedures_main_scripts/general_scripts/second_step_lambda_adapt_function_lsodar_of_grhmc_adaptive.R")
source("tuning_procedures_main_scripts/general_scripts/third_step_general_function_of_grhmc_adaptive.R")

final_VARI_general_function_of_grhmc_adaptive <- function(
    log_target_grad,
    lambda_initial,
    n_parameters = NULL,
    time_period_adaptive = 5000,
    time_period_generating_samples = 5000,
    n_generated_samples = 5000,
    n_adaptive_samples = 1000,
    diag_s_elements_initial = NULL,
    m_initial = NULL,
    qbar_initial = NULL,
    pbar_initial = NULL,
    random_state = NULL,
    maxsteps = NULL,
    rtol = NULL,
    atol = NULL,
    proportion_time_until_adaptive_start = 0.05,
    min_proportion_of_previous_state = 0.5, # Only applies for VARI or ISG
    max_proportion_of_previous_state = 2,
    lambda_adaptive_step_sim_time = 5000,
    lambda_adaptive_step_ema_beta = 0.99,
    lambda_adaptive_step_factor = 1,
    lambda_adaptive_step_n_samples = 1000,
    lambda_adaptive_step_max_nut_time = 100
){
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  first_step_run <- first_step_grhmc_transformed_function_VARI(
    log_target_grad = log_target_grad,
    lambda = lambda_initial,
    n_parameters = n_parameters,
    T = time_period_adaptive,
    n_samples = n_adaptive_samples,
    diag_s_elements_initial = diag_s_elements_initial,
    m_initial = m_initial,
    qbar_initial = qbar_initial,
    pbar_initial = pbar_initial,
    maxsteps = maxsteps,
    rtol = rtol,
    atol = atol,
    proportion_time_until_adaptive_start = proportion_time_until_adaptive_start,
    min_proportion_of_previous_state = min_proportion_of_previous_state,
    max_proportion_of_previous_state = max_proportion_of_previous_state
  )
  
  second_step_run <- second_step_lambda_adapt_function(
    fun = log_target_grad,
    q0 = first_step_run$qbar_end,
    p0 = first_step_run$pbar_end,
    m.c = first_step_run$m_adaptive,
    S.c = first_step_run$s_adaptive,
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
    T = time_period_generating_samples,
    n_samples = n_generated_samples,
    diag_s_elements_initial = first_step_run$s_adaptive,
    m_initial = first_step_run$m_adaptive,
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

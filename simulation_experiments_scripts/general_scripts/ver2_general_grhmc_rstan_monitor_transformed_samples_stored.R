# RSTAN GENERAL MONITOR EXAMPLE
rm(list = ls())
# Load the general functions to simulate
source("tuning_procedures_main_scripts/VARI/final_VARI_general_function_of_grhmc_adaptive.R")
source("tuning_procedures_main_scripts/MCT/ver2_final_MCT_general_function_of_grhmc_adaptive.R")
source("tuning_procedures_main_scripts/ISG/final_ISG_general_function_of_grhmc_adaptive.R")

general_grhmc_rstan_monitor_samples_stored <- function(
  export_result = "grhmc_monitor_output.RDS",
  log_target_grad,
  lambda_initial = 0.2,
  n_chains = 10,
  n_generated_samples = 1000,
  time_period_generating_samples = 5000,
  n_parameters,
  desired_time_until_crossing_median = pi,
  diag_s_elements_initial = NULL,
  m_initial = NULL,
  qbar_initial = NULL,
  pbar_initial = NULL,
  maxsteps = 5000,
  rtol = 1e-6,
  atol = 1e-6,
  random_state = NULL,
  vari_and_isg_time_period_adaptive = 5000,
  vari_and_isg_n_adaptive_samples = 1000,
  vari_and_isg_proportion_time_until_adaptive_start = 0.05,
  vari_and_isg_min_proportion_of_previous_state = 0.5,
  vari_and_isg_max_proportion_of_previous_state = 2,
  mct_step_sim_time_first_run = 1000,
  mct_step_sim_time_second_run = 5000,
  mct_step_proportion_time_until_adaptive = 0.05,
  mct_step_time_step_update_data_median = 1,
  mct_step_n_adaptive_samples = 1000,
  mct_step_kappa_elements_dual_averaging = NULL,
  mct_step_gamma_elements_dual_averaging_first_run = NULL,
  mct_step_gamma_elements_dual_averaging_second_run = NULL,
  mct_step_t0_elements_dual_averaging = NULL,
  mct_step_mu_elements_dual_averaging_first_run = NULL,
  mct_step_mu_log_scale_second_run = TRUE,
  mct_step_mu_scale_factor_second_run = 1.1,
  lambda_adaptive_step_sim_time = 5000,
  lambda_adaptive_step_ema_beta = 0.99,
  lambda_adaptive_step_factor = 1,
  lambda_adaptive_step_n_samples = 1000,
  lambda_adaptive_step_max_nut_time = 100
){
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }

  mct_step_time_period_adaptive <- mct_step_sim_time_first_run + mct_step_sim_time_second_run
  time_period_adaptive <- mct_step_time_period_adaptive
  
  if(abs(mct_step_time_period_adaptive - vari_and_isg_time_period_adaptive) >= 1e-7){
    stop("Make sure that the length of the adaptive time period is the same between the methods for comparability!")
  }
  
  # Nr.1 - Run 10 iterations to get 10 chains using the ISG
  
  transformed_ISG_3d_array <- array(
    # first dimension: Number of samples
    # second dimension: Number of trajectories/chains 
    # third dimension: Number of parameters/dimension of position variable
    dim = c(n_generated_samples, n_chains, n_parameters) 
  )
  
  transformed_ISG_n_evals_ode <- rep(0, n_chains)
  transformed_ISG_lambda_adaptive <- rep(0, n_chains)
  
  transformed_ISG_m_elements <- matrix(nrow = n_chains, ncol = n_parameters)
  transformed_ISG_s_elements <- matrix(nrow = n_chains, ncol = n_parameters)
  
  for(i in 1:n_chains){
    
    transformed_ISG_run <- final_ISG_general_function_of_grhmc_adaptive(
      
      log_target_grad = log_target_grad,
      lambda_initial = lambda_initial,
      n_parameters = n_parameters,
      time_period_adaptive = vari_and_isg_time_period_adaptive,
      time_period_generating_samples = time_period_generating_samples,
      n_generated_samples = n_generated_samples,
      n_adaptive_samples = vari_and_isg_n_adaptive_samples,
      diag_s_elements_initial = diag_s_elements_initial,
      m_initial = m_initial,
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      maxsteps = maxsteps,
      rtol = rtol,
      atol = atol,
      proportion_time_until_adaptive_start = vari_and_isg_proportion_time_until_adaptive_start,
      min_proportion_of_previous_state = vari_and_isg_min_proportion_of_previous_state,
      max_proportion_of_previous_state = vari_and_isg_max_proportion_of_previous_state,
      lambda_adaptive_step_sim_time = lambda_adaptive_step_sim_time,
      lambda_adaptive_step_ema_beta = lambda_adaptive_step_ema_beta,
      lambda_adaptive_step_factor = lambda_adaptive_step_factor,
      lambda_adaptive_step_n_samples = lambda_adaptive_step_n_samples,
      lambda_adaptive_step_max_nut_time = lambda_adaptive_step_max_nut_time 
    )
    
    name <- paste("transformed_ISG_run_nr", sep = "_", i)
    assign(name, transformed_ISG_run$final_step_run$q_original_samples)
    
    transformed_ISG_3d_array[, i, ] <- eval(parse(text = name))
    
    transformed_ISG_n_evals_ode[i] <- transformed_ISG_run$final_step_run$n_evals_ode
    transformed_ISG_lambda_adaptive[i] <- transformed_ISG_run$final_step_run$lambda
    
    transformed_ISG_m_elements[i, ] <- transformed_ISG_run$final_step_run$m_elements
    transformed_ISG_s_elements[i, ] <- transformed_ISG_run$final_step_run$s_elements 
    
    print(paste0("ISG nr. ", i))
    
  }
  
  # Check that each "slice" (i.e. third dimension) corresponds to a given position variable
  
  # colMeans(transformed_ISG_3d_array[, , 1]) 
  # apply(transformed_ISG_3d_array[, , 1], MARGIN = 2, sd) 
  # apply(transformed_ISG_3d_array[, , 1], MARGIN = 2, var) 
  # 
  # colMeans(transformed_ISG_3d_array[, , 2]) 
  # apply(transformed_ISG_3d_array[, , 2], MARGIN = 2, sd) 
  # apply(transformed_ISG_3d_array[, , 2], MARGIN = 2, var) 
  
  # Run this 3d-array into rstan::monitor, all these samples are samples after burn in period
  
  transformed_ISG_mc_summary <- rstan::monitor(
    
    sims = transformed_ISG_3d_array,
    warmup = 0
    
  )
  
  transformed_ISG_mc_summary
  
  # Nr.2 - Now do the same using VARI
  
  transformed_VARI_3d_array <- array(
    # first dimension: Number of samples
    # second dimension: Number of trajectories/chains 
    # third dimension: Number of parameters/dimension of position variable
    dim = c(n_generated_samples, n_chains, n_parameters) 
  )
  
  transformed_VARI_n_evals_ode <- rep(0, n_chains)
  transformed_VARI_lambda_adaptive <- rep(0, n_chains)
  
  transformed_VARI_m_elements <- matrix(nrow = n_chains, ncol = n_parameters)
  transformed_VARI_s_elements <- matrix(nrow = n_chains, ncol = n_parameters)
  
  for(i in 1:n_chains){
    
    transformed_VARI_run <- final_VARI_general_function_of_grhmc_adaptive(
      
      log_target_grad = log_target_grad,
      lambda_initial = lambda_initial,
      n_parameters = n_parameters,
      time_period_adaptive = vari_and_isg_time_period_adaptive,
      time_period_generating_samples = time_period_generating_samples,
      n_generated_samples = n_generated_samples,
      n_adaptive_samples = vari_and_isg_n_adaptive_samples,
      diag_s_elements_initial = diag_s_elements_initial,
      m_initial = m_initial,
      qbar_initial = qbar_initial,
      pbar_initial = pbar_initial,
      maxsteps = maxsteps,
      rtol = rtol,
      atol = atol,
      proportion_time_until_adaptive_start = vari_and_isg_proportion_time_until_adaptive_start,
      min_proportion_of_previous_state = vari_and_isg_min_proportion_of_previous_state,
      max_proportion_of_previous_state = vari_and_isg_max_proportion_of_previous_state,
      lambda_adaptive_step_sim_time = lambda_adaptive_step_sim_time,
      lambda_adaptive_step_ema_beta = lambda_adaptive_step_ema_beta,
      lambda_adaptive_step_factor = lambda_adaptive_step_factor,
      lambda_adaptive_step_n_samples = lambda_adaptive_step_n_samples,
      lambda_adaptive_step_max_nut_time = lambda_adaptive_step_max_nut_time 
    )
    
    name <- paste("transformed_VARI_run_nr", sep = "_", i)
    assign(name, transformed_VARI_run$final_step_run$q_original_samples)
    
    transformed_VARI_3d_array[, i, ] <- eval(parse(text = name))
    
    transformed_VARI_n_evals_ode[i] <- transformed_VARI_run$final_step_run$n_evals_ode
    transformed_VARI_lambda_adaptive[i] <- transformed_VARI_run$final_step_run$lambda
      
    transformed_VARI_m_elements[i, ] <- transformed_VARI_run$final_step_run$m_elements
    transformed_VARI_s_elements[i, ] <- transformed_VARI_run$final_step_run$s_elements 
    
    print(paste0("VARI nr. ", i))
  }
  
  # Check that each "slice" (i.e. third dimension) corresponds to a given position variable
  
  # colMeans(transformed_VARI_3d_array[, , 1]) 
  # apply(transformed_VARI_3d_array[, , 1], MARGIN = 2, sd) 
  # apply(transformed_VARI_3d_array[, , 1], MARGIN = 2, var) 
  # 
  # 
  # colMeans(transformed_VARI_3d_array[, , 2]) 
  # apply(transformed_VARI_3d_array[, , 2], MARGIN = 2, sd) 
  # apply(transformed_VARI_3d_array[, , 2], MARGIN = 2, var) 
  
  
  # Run this 3d-array into rstan::monitor, all these samples are samples after burn in period
  
  transformed_VARI_mc_summary <- rstan::monitor(
    
    sims = transformed_VARI_3d_array,
    warmup = 0
    
  )
  
  transformed_VARI_mc_summary
  
  # Nr.3 - Now do the same using median crossing times
  
  transformed_MCT_3d_array <- array(
    # first dimension: Number of samples
    # second dimension: Number of trajectories/chains 
    # third dimension: Number of parameters/dimension of position variable
    dim = c(n_generated_samples, n_chains, n_parameters) 
  )
  
  transformed_MCT_n_evals_ode <- rep(0, n_chains)
  transformed_MCT_lambda_adaptive <- rep(0, n_chains)
  
  transformed_MCT_m_elements <- matrix(nrow = n_chains, ncol = n_parameters)
  transformed_MCT_s_elements <- matrix(nrow = n_chains, ncol = n_parameters)
  
  for(i in 1:n_chains){
    
    transformed_MCT_run <- final_MCT_general_function_of_grhmc_adaptive(
      
      log_target_grad = log_target_grad,
      lambda_initial = lambda_initial,
      n_parameters = n_parameters,
      median_crossing_step_sim_time_first_run = mct_step_sim_time_first_run,
      median_crossing_step_sim_time_second_run = mct_step_sim_time_second_run,
      median_crossing_step_n_samples = mct_step_n_adaptive_samples,
      desired_time_until_crossing_median = desired_time_until_crossing_median,
      median_crossing_step_diag_s_elements_initial = diag_s_elements_initial,
      median_crossing_step_m_initial = m_initial,
      median_crossing_step_qbar_initial = qbar_initial,
      median_crossing_step_pbar_initial = pbar_initial,
      maxsteps = maxsteps,
      rtol = rtol,
      atol = atol,
      median_crossing_step_proportion_time_until_adaptive = mct_step_proportion_time_until_adaptive,
      median_crossing_step_time_step_update_data_median = mct_step_time_step_update_data_median,
      median_crossing_step_kappa_elements_dual_averaging = mct_step_kappa_elements_dual_averaging,
      median_crossing_step_gamma_elements_dual_averaging_first_run = mct_step_gamma_elements_dual_averaging_first_run,
      median_crossing_step_gamma_elements_dual_averaging_second_run = mct_step_gamma_elements_dual_averaging_second_run,
      median_crossing_step_t0_elements_dual_averaging = mct_step_t0_elements_dual_averaging,
      median_crossing_step_mu_elements_dual_averaging_first_run = mct_step_mu_elements_dual_averaging_first_run,
      median_crossing_step_mu_log_scale_second_run = mct_step_mu_log_scale_second_run,
      median_crossing_step_mu_scale_factor_second_run = mct_step_mu_scale_factor_second_run,
      lambda_adaptive_step_sim_time = lambda_adaptive_step_sim_time,
      lambda_adaptive_step_ema_beta = lambda_adaptive_step_ema_beta,
      lambda_adaptive_step_factor = lambda_adaptive_step_factor,
      lambda_adaptive_step_n_samples = lambda_adaptive_step_n_samples,
      lambda_adaptive_step_max_nut_time = lambda_adaptive_step_max_nut_time,
      final_sampling_step_sim_time = time_period_generating_samples,
      final_sampling_step_n_samples = n_generated_samples
      
    )
    
    name <- paste("transformed_MCT_run_nr", sep = "_", i)
    assign(name, transformed_MCT_run$final_step_run$q_original_samples)
    
    transformed_MCT_3d_array[, i, ] <- eval(parse(text = name))
    
    transformed_MCT_n_evals_ode[i] <- transformed_MCT_run$final_step_run$n_evals_ode
    transformed_MCT_lambda_adaptive[i] <- transformed_MCT_run$final_step_run$lambda
    
    transformed_MCT_m_elements[i, ] <- transformed_MCT_run$final_step_run$m_elements
    transformed_MCT_s_elements[i, ] <- transformed_MCT_run$final_step_run$s_elements 
    
    print(paste0("MCT nr. ", i))
    
  }
  
  # Check that each "slice" (i.e. third dimension) corresponds to a given position variable
  
  # colMeans(transformed_MCT_3d_array[, , 1]) 
  # apply(transformed_MCT_3d_array[, , 1], MARGIN = 2, sd) 
  # apply(transformed_MCT_3d_array[, , 1], MARGIN = 2, var) 
  # 
  # 
  # colMeans(transformed_MCT_3d_array[, , 2]) 
  # apply(transformed_MCT_3d_array[, , 2], MARGIN = 2, sd) 
  # apply(transformed_MCT_3d_array[, , 2], MARGIN = 2, var) 
  
  
  # Run this 3d-array into rstan::monitor, all these samples are samples after burn in period
  
  transformed_MCT_mc_summary <- rstan::monitor(
    
    sims = transformed_MCT_3d_array,
    warmup = 0
    
  )
  
  transformed_MCT_mc_summary
  
  # Comparison
  
  summary_VARI <- transformed_VARI_mc_summary
  summary_ISG <- transformed_ISG_mc_summary
  summary_MCT <- transformed_MCT_mc_summary
  
  n_eff_VARI <- transformed_VARI_mc_summary$n_eff
  n_eff_ISG <- transformed_ISG_mc_summary$n_eff
  n_eff_MCT <- transformed_MCT_mc_summary$n_eff
  
  n_evals_ode_VARI <- sum(transformed_VARI_n_evals_ode)
  n_evals_ode_ISG <- sum(transformed_ISG_n_evals_ode)
  n_evals_ode_MCT <- sum(transformed_MCT_n_evals_ode)
  
  n_eff_per_eval_ode_VARI <- n_eff_VARI / n_evals_ode_VARI
  n_eff_per_eval_ode_ISG <- n_eff_ISG / n_evals_ode_ISG
  n_eff_per_eval_ode_MCT <- n_eff_MCT / n_evals_ode_MCT
  
  summary_table <- data.frame(
    method = c(rep("VARI", n_parameters), rep("ISG", n_parameters), rep("Median crossing", n_parameters)),
    position_coordinate = rep(1:n_parameters, 3),
    n_eff = c(n_eff_VARI, n_eff_ISG, n_eff_MCT),
    n_evals_ode = c(rep(n_evals_ode_VARI, n_parameters), rep(n_evals_ode_ISG, n_parameters), rep(n_evals_ode_MCT, n_parameters)),
    n_eff_per_eval_ode = c(n_eff_per_eval_ode_VARI, n_eff_per_eval_ode_ISG, n_eff_per_eval_ode_MCT)
  ) 
  
  return_list <- list(
    time_period = time_period_adaptive + lambda_adaptive_step_sim_time + time_period_generating_samples,
    time_period_adaptive = time_period_adaptive,
    lambda_adaptive_step_sim_time = lambda_adaptive_step_sim_time,
    time_period_generating_samples = time_period_generating_samples,
    n_generated_samples = n_generated_samples,
    n_chains = n_chains,
    VARI_3d_array = transformed_VARI_3d_array,
    ISG_3d_array = transformed_ISG_3d_array,
    MCT_3d_array = transformed_MCT_3d_array,
    summary_VARI = summary_VARI,
    summary_ISG = summary_ISG,
    summary_MCT = summary_MCT,
    n_eff_VARI = n_eff_VARI,
    n_eff_ISG = n_eff_ISG,
    n_eff_MCT = n_eff_MCT,
    n_evals_ode_VARI = n_evals_ode_VARI,
    n_evals_ode_ISG = n_evals_ode_ISG,
    n_evals_ode_MCT = n_evals_ode_MCT,
    n_eff_per_eval_ode_VARI = n_eff_per_eval_ode_VARI,
    n_eff_per_eval_ode_ISG = n_eff_per_eval_ode_ISG,
    n_eff_per_eval_ode_MCT = n_eff_per_eval_ode_MCT,
    VARI_s_elements = transformed_VARI_s_elements,
    ISG_s_elements = transformed_ISG_s_elements,
    MCT_s_elements = transformed_MCT_s_elements,
    VARI_m_elements = transformed_VARI_m_elements,
    ISG_m_elements = transformed_ISG_m_elements,
    MCT_m_elements = transformed_MCT_m_elements,
    VARI_lambda_adaptive = transformed_VARI_lambda_adaptive,
    ISG_lambda_adaptive = transformed_ISG_lambda_adaptive,
    MCT_lambda_adaptive = transformed_MCT_lambda_adaptive,
    summary_table = summary_table,
    maxsteps = maxsteps,
    rtol = rtol,
    atol = atol,
    random_state = random_state
  )
  
  saveRDS(return_list, export_result, compress = FALSE)
  
  return(return_list)
  
}


# First step - Tuning m and S using MCT

# Split in two runs as mentioned in Appendix

source("tuning_procedures_main_scripts/MCT/main_procedure_MCT_general_function_of_grhmc_adaptive.R")

first_step_grhmc_transformed_function_MCT <- function(log_target_grad,
                                                                        lambda,
                                                                        n_parameters = NULL,
                                                                        sim_time_first_run = 1000,
                                                                        sim_time_second_run = 5000,
                                                                        n_samples = 1000,
                                                                        desired_time_until_crossing_median = pi,
                                                                        diag_s_elements_initial = NULL, 
                                                                        m_initial = NULL,
                                                                        qbar_initial = NULL,
                                                                        pbar_initial = NULL,
                                                                        random_state = NULL,
                                                                        maxsteps = NULL,
                                                                        rtol = NULL,
                                                                        atol = NULL,
                                                                        proportion_time_until_adaptive_start = 0.05,
                                                                        time_step_update_data_median = 1,
                                                                        kappa_elements_dual_averaging = NULL,
                                                                        gamma_elements_dual_averaging_first_run = NULL,
                                                                        gamma_elements_dual_averaging_second_run = NULL,
                                                                        t0_elements_dual_averaging = NULL,
                                                                        mu_elements_dual_averaging_first_run = NULL,
                                                                        mu_log_scale_second_run = TRUE,
                                                                        mu_scale_factor_second_run = NULL
){
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  n_dim_error <- FALSE
  
  tryCatch(
    {
      n_dim <- length(log_target_grad(0))
    },
    error = function(e){
      print(e)
      n_dim_error <<- TRUE
    }
  )
  
  if(n_dim_error == TRUE){
    if(is.null(n_parameters)){
      stop("Please specify n_parameters")
    } else {
      n_dim <- n_parameters
    }
  }  
  
  if(is.null(gamma_elements_dual_averaging_first_run)){
    gamma_elements_dual_averaging_first_run <- rep(10, n_dim)
  }
  
  if(is.null(gamma_elements_dual_averaging_second_run)){
    gamma_elements_dual_averaging_second_run <- rep(25, n_dim)
  }
  
  if(is.null(mu_scale_factor_second_run)){
    mu_scale_factor_second_run <- rep(1.1, n_dim)
  }
  
  first_run <- main_procedure_grhmc_transformed_function_MCT(
    log_target_grad = log_target_grad,
    lambda = lambda,
    n_parameters = n_parameters,
    T = sim_time_first_run,
    n_samples = n_samples,
    desired_time_until_crossing_median = desired_time_until_crossing_median,
    diag_s_elements_initial = diag_s_elements_initial,
    m_initial = m_initial,
    qbar_initial = qbar_initial,
    pbar_initial = pbar_initial,
    maxsteps = maxsteps,
    rtol = rtol,
    atol = atol,
    proportion_time_until_adaptive_start = proportion_time_until_adaptive_start,
    time_step_update_data_median = time_step_update_data_median,
    kappa_elements_dual_averaging = kappa_elements_dual_averaging,
    gamma_elements_dual_averaging = gamma_elements_dual_averaging_first_run,
    t0_elements_dual_averaging = t0_elements_dual_averaging,
    mu_elements_dual_averaging = mu_elements_dual_averaging_first_run
  )
  
  s_elements_from_first_run <- first_run$s_adaptive
  m_elements_from_first_run <- first_run$m_adaptive
  qbar_end_from_first_run <- first_run$qbar_end
  pbar_end_from_first_run <- first_run$pbar_end
  
  if(mu_log_scale_second_run){
    mu_elements_dual_averaging_second_run <- mu_scale_factor_second_run * log(s_elements_from_first_run)
  } else {
    mu_elements_dual_averaging_second_run <- log(mu_scale_factor_second_run * s_elements_from_first_run)
  }
  
  second_run <- main_procedure_grhmc_transformed_function_MCT(
    log_target_grad = log_target_grad,
    lambda = lambda,
    n_parameters = n_parameters,
    T = sim_time_second_run,
    n_samples = n_samples,
    desired_time_until_crossing_median = desired_time_until_crossing_median,
    diag_s_elements_initial = s_elements_from_first_run,
    m_initial = m_elements_from_first_run,
    qbar_initial = qbar_end_from_first_run,
    pbar_initial = pbar_end_from_first_run,
    maxsteps = maxsteps,
    rtol = rtol,
    atol = atol,
    proportion_time_until_adaptive_start = proportion_time_until_adaptive_start,
    time_step_update_data_median = time_step_update_data_median,
    kappa_elements_dual_averaging = kappa_elements_dual_averaging,
    gamma_elements_dual_averaging = gamma_elements_dual_averaging_second_run,
    t0_elements_dual_averaging = t0_elements_dual_averaging,
    mu_elements_dual_averaging = mu_elements_dual_averaging_second_run
  )
  
  return(
    list(
      first_median_crossing_run = first_run,
      second_median_crossing_run = second_run
    )
  )
  
}

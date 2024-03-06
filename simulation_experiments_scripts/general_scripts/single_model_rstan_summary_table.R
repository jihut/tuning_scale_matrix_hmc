rstan_summary_table <- function(RDS_object, sample_grid, name, round_decimals_n_eff_per_100000_eval_ode = 0, round_decimals_mean_s_elements = 2, round_decimals_sd_s_elements = 3){
  
  dim_3d_array <- dim(RDS_object$VARI_3d_array)
  
  n_parameters <- dim_3d_array[3]
  n_samples <- dim_3d_array[1]
  
  if((n_samples %% sample_grid) != 0){
    stop("Need to choose sample_grid that divides the number of samples per chain!")
  }
  
  transformed_VARI_mc_summary <- rstan::monitor(
    sims = RDS_object$VARI_3d_array[seq(from = sample_grid, to = n_samples, by = sample_grid), , ],
    warmup = 0
  )
  
  transformed_ISG_mc_summary <- rstan::monitor(
    sims = RDS_object$ISG_3d_array[seq(from = sample_grid, to = n_samples, by = sample_grid), , ],
    warmup = 0
  )
  
  transformed_MCT_mc_summary <- rstan::monitor(
    sims = RDS_object$MCT_3d_array[seq(from = sample_grid, to = n_samples, by = sample_grid), , ],
    warmup = 0
  )
  
  n_eff_VARI <- transformed_VARI_mc_summary$n_eff
  n_eff_ISG <- transformed_ISG_mc_summary$n_eff
  n_eff_MCT <- transformed_MCT_mc_summary$n_eff
  
  n_evals_ode_VARI <- RDS_object$n_evals_ode_VARI
  n_evals_ode_ISG <- RDS_object$n_evals_ode_ISG
  n_evals_ode_MCT <- RDS_object$n_evals_ode_MCT
  
  n_eff_per_eval_ode_VARI <- n_eff_VARI / n_evals_ode_VARI
  n_eff_per_eval_ode_ISG <- n_eff_ISG / n_evals_ode_ISG
  n_eff_per_eval_ode_MCT <- n_eff_MCT / n_evals_ode_MCT
  
  min_max_indices_VARI <- c(which.min(n_eff_VARI), which.max(n_eff_VARI))
  min_max_indices_ISG <- c(which.min(n_eff_ISG), which.max(n_eff_ISG))
  min_max_indices_MCT <- c(which.min(n_eff_MCT), which.max(n_eff_MCT))
  
  mean_s_elements_VARI <- format(round(apply(RDS_object$VARI_s_elements, 2, mean), round_decimals_mean_s_elements), nsmall = round_decimals_mean_s_elements)
  sd_s_elements_VARI <- format(round(apply(RDS_object$VARI_s_elements, 2, sd), round_decimals_sd_s_elements), nsmall = round_decimals_sd_s_elements)
  mean_sd_s_elements_VARI <- paste(mean_s_elements_VARI, " (", sd_s_elements_VARI, ")", sep = "")
  
  mean_s_elements_ISG <- format(round(apply(RDS_object$ISG_s_elements, 2, mean), round_decimals_mean_s_elements), nsmall = round_decimals_mean_s_elements)
  sd_s_elements_ISG <- format(round(apply(RDS_object$ISG_s_elements, 2, sd), round_decimals_sd_s_elements), nsmall = round_decimals_sd_s_elements)
  mean_sd_s_elements_ISG <- paste(mean_s_elements_ISG, " (", sd_s_elements_ISG, ")", sep = "")
  
  mean_s_elements_MCT <- format(round(apply(RDS_object$MCT_s_elements, 2, mean), round_decimals_mean_s_elements), nsmall = round_decimals_mean_s_elements)
  sd_s_elements_MCT <- format(round(apply(RDS_object$MCT_s_elements, 2, sd), round_decimals_sd_s_elements), nsmall = round_decimals_sd_s_elements)
  mean_sd_s_elements_MCT <- paste(mean_s_elements_MCT, " (", sd_s_elements_MCT, ")", sep = "")
  
  summary_table <- data.frame(
    name = name,
    method = c(rep("VARI", n_parameters), rep("ISG", n_parameters), rep("MCT", n_parameters)),
    position_coordinate = rep(1:n_parameters, 3),
    mean_sd_s_elements = c(mean_sd_s_elements_VARI, mean_sd_s_elements_ISG, mean_sd_s_elements_MCT),
    n_eff = c(n_eff_VARI, n_eff_ISG, n_eff_MCT),
    n_evals_ode = c(rep(n_evals_ode_VARI, n_parameters), rep(n_evals_ode_ISG, n_parameters), rep(n_evals_ode_MCT, n_parameters)),
    n_eff_per_100000_eval_ode = round(c(n_eff_per_eval_ode_VARI, n_eff_per_eval_ode_ISG, n_eff_per_eval_ode_MCT) * 100000, round_decimals_n_eff_per_100000_eval_ode)
  ) 
  
  min_max_summary_table <- data.frame(
    name = name,
    method = c(rep("VARI", 2), rep("ISG", 2), rep("MCT", 2)),
    position_coordinate = rep(c("min ESS", "max ESS"), 3),
    mean_sd_s_elements = c(mean_sd_s_elements_VARI[min_max_indices_VARI], mean_sd_s_elements_ISG[min_max_indices_ISG], mean_sd_s_elements_MCT[min_max_indices_MCT]),
    n_eff = c(n_eff_VARI[min_max_indices_VARI], n_eff_ISG[min_max_indices_ISG], n_eff_MCT[min_max_indices_MCT]),
    n_evals_ode = c(rep(n_evals_ode_VARI, 2), rep(n_evals_ode_ISG, 2), rep(n_evals_ode_MCT, 2)),
    n_eff_per_100000_eval_ode = round(c(n_eff_per_eval_ode_VARI[min_max_indices_VARI], n_eff_per_eval_ode_ISG[min_max_indices_ISG], n_eff_per_eval_ode_MCT[min_max_indices_MCT]) * 100000, round_decimals_n_eff_per_100000_eval_ode)
  ) 
  
  list(all_summary_table = summary_table, min_max_summary_table = min_max_summary_table)
  
}

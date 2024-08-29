# Main function for tuning m and S using median crossing times with GRHMC framework

source("tuning_procedures_main_scripts/general_scripts/p2_algorithm_class.R")
source("tuning_procedures_main_scripts/general_scripts/dual_averaging_class.R")

main_procedure_grhmc_transformed_function_MCT <- function(log_target_grad,
                                                          lambda,
                                                          n_parameters = NULL,
                                                          T = 5000,
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
                                                          gamma_elements_dual_averaging = NULL,
                                                          t0_elements_dual_averaging = NULL,
                                                          mu_elements_dual_averaging = NULL
){
  
  if(proportion_time_until_adaptive_start * T / time_step_update_data_median < 5){
    stop("The chosen time step to update median and the proportion of time until the start adaptive does not generate 5 samples before the start adaptive time point! 
         Pick a smaller time step or increase the proportion!")
  }
  
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
  
  n_evals_ode <- 0L # store how many times the ode function below is called in the process
  
  count_n_evals_ode <- function(){
    
    n_evals_ode <<- n_evals_ode + 1 # increment the counter every time the ode function is called
    
  }
  
  if(is.null(m_initial)){
    m_adaptive <- rep(0, n_dim) # initialize m 
  } else {
    m_adaptive <- m_initial
  }
  old_m_adaptive <- m_adaptive
  
  if(is.null(diag_s_elements_initial)){
    s_adaptive <- rep(1, n_dim) # initialize m 
  } else {
    s_adaptive <- diag_s_elements_initial
  }
  old_s_adaptive <- s_adaptive
  
  if(is.null(kappa_elements_dual_averaging)){
    kappa_elements_dual_averaging <- rep(0.75, n_dim)
  } else {
    kappa_elements_dual_averaging <- kappa_elements_dual_averaging
  }
  
  if(is.null(gamma_elements_dual_averaging)){
    # gamma_elements_dual_averaging <- rep(0.5, n_dim)
    gamma_elements_dual_averaging <- rep(10, n_dim)
  } else {
    gamma_elements_dual_averaging <- gamma_elements_dual_averaging
  }
  
  if(is.null(t0_elements_dual_averaging)){
    t0_elements_dual_averaging <- rep(10, n_dim)
    # t0_elements_dual_averaging <- rep(5, n_dim)
  } else {
    t0_elements_dual_averaging <- t0_elements_dual_averaging
  }
  
  if(is.null(mu_elements_dual_averaging)){
    mu_elements_dual_averaging <- rep(0, n_dim)
  } else {
    mu_elements_dual_averaging <- mu_elements_dual_averaging
  }
  
  n_q_root_median <- rep(0, n_dim) # store number of times a given position coordinate has crossed median value
  
  # diff_time_q_root_median <- rep_len(list(numeric()), n_dim) # store the time between two events of crossing median for each position coordinate
  
  # time_q_root_median <- rep_len(list(numeric()), n_dim) # store the specific time points where the median of a given position coordinate has been crossed
  
  cross_median_time <- rep(0, n_dim) # store the fictitious time point when the median is crossed, this will be used to find the time between such events
  
  if(is.null(qbar_initial)){
    # qbar_initial <- rep(0, n_dim)
    qbar_initial <- rnorm(n_dim) # If initial values of qbar are not specified --> simulate randomly from a multivariate normal distribution with zero mean vector and diagonal covariance matrix
  } else {
    qbar_initial <- qbar_initial
  }
  
  
  if(!is.null(maxsteps)){
    maxsteps <- maxsteps
  } else {
    maxsteps <- 5000 # default in deSolve::lsodar
  }
  
  if(is.null(rtol)){
    rtol <- 1e-6 # default in deSolve::lsodar
  } else {
    rtol <- rtol
  }
  
  if(is.null(atol)){
    atol <- 1e-6 # default in deSolve::lsodar
  } else {
    atol <- atol
  }
  
  
  # define the (first order) ode of state
  
  # state[1]: number of events up to time point
  # state[2]: Lambda, which resets after each event
  # state[3]: u, simulate from exp(1) after each event
  # state[4]: time to update data for median estimation
  # state[5:(5 + n_dim - 1)]: qbar 
  # state[(5 + n_dim):(5 + 2*n_dim - 1)]: pbar
  # state[(5 + 2*n_dim):(5 + 3*n_dim - 1)]: \int qbar dt
  # state[(5 + 3*n_dim):(5 + 4*n_dim - 1)]: \int qbar^2 dt
  # state[(5 + 4*n_dim):(5 + 5*n_dim - 1)]: m - Median adaptive - stops changing after T/2 time units
  # state[(5 + 5*n_dim):(5 + 6*n_dim - 1)]: S - diagonal elements of S adaptive - stops changing after T/2 time units
  
  # Initialize P2 algorithm class for each coordinate as we need the median for each position coordinate
  
  for(i in 1:n_dim){
    
    assign(paste("q", i, "_median", sep = ""), p2_algorithm_class$new())
    assign(
      paste("s", i, "_solution", sep = ""), 
      dual_averaging_class$new(
        kappa = kappa_elements_dual_averaging[i],
        gamma = gamma_elements_dual_averaging[i],
        t0 = t0_elements_dual_averaging[i],
        mu = mu_elements_dual_averaging[i],
        x_initial = log(s_adaptive[i])
      )
    ) # notice log scale here
    
  }
  
  # Define ODE
  
  ode <- function(t, state, parms){
    
    count_n_evals_ode()
    
    ret <- c(
      
      0, # number of events
      lambda, 
      0, # u
      0, # time to update data for median estimation - only change after an update
      state[(5 + n_dim):(5 + 2*n_dim - 1)], # dot qbar = pbar
      state[(5 + 5*n_dim):(5 + 6*n_dim - 1)] * 
        log_target_grad(state[(5 + 4*n_dim):(5 + 5*n_dim - 1)] + state[(5 + 5*n_dim):(5 + 6*n_dim - 1)] * state[5:(5 + n_dim - 1)]), # dot pbar = gradient of log transformed density wrt. qbar
      state[5:(5 + n_dim - 1)], # int qbar dt
      state[5:(5 + n_dim - 1)]^2, # int qbar^2 dt
      # state[(5 + 4*n_dim):(5 + 5*n_dim - 1)] + state[(5 + 5*n_dim):(5 + 6*n_dim - 1)] * state[5:(5 + n_dim - 1)], # int q dt
      # (state[(5 + 4*n_dim):(5 + 5*n_dim - 1)] + state[(5 + 5*n_dim):(5 + 6*n_dim - 1)] * state[5:(5 + n_dim - 1)])^2, # int q^2 dt
      rep(0, n_dim), # median adaptive - only update at events
      rep(0, n_dim) # diagonal elements of S adaptive - only update at events
    )
    
    return(list(ret))
    
  }
  
  root.fun <- function(t, state, parms){
    
    return(
      
      c(
        state[2] - state[3], # root that leads to standard momentum refresh
        t - state[4], # root that leads to adding a new sample to the sample vector used for median estimation
        t - (proportion_time_until_adaptive_start * T - 1e-9), # add a single event at this time point so that median can be updated at this point in case the position crosses the initial median before an event has happened - t
        state[5:(5 + n_dim - 1)] # events whenever a coordinate crosses its corresponding median value --> qbar crosses 0
      )
      
    )
    
  }
  
  at.event <- function(t, state, parms){
    
    yroot <- c(
      state[2] - state[3], # root that leads to standard momentum refresh
      t - state[4], # root that leads to adding a new sample to the sample vector used for median estimation
      t - (proportion_time_until_adaptive_start * T - 1e-9), # add a single event at this time point so that median can be updated at this point in case the position crosses the initial median before an event has happened - t
      state[5:(5 + n_dim - 1)] # events whenever a coordinate crosses its corresponding median value --> qbar crosses 0
    )
    
    whichroot <- which(min(abs(yroot)) == abs(yroot))
    
    if(t > 0){ # when qbar_initial is equal to initialized median, we get a root right at the start. To avoid this, we set this condition here. 
      
      if(whichroot == 1){ # standard momentum refresh event
        
        
        p_refresh <- rnorm(n_dim)
        u_refresh <- rexp(1)
        
        newState <- c(
          
          state[1] + 1, # update the number of momentum refresh event counter
          0, # reset Lambda
          u_refresh, # simulate a new u from exponential distribution with rate 1
          state[4], # time to update data for median estimation
          state[5:(5 + n_dim - 1)], # qbar unchanged since m and S are now unchanged during a momentum refresh event
          p_refresh, # refresh momentum
          state[(5 + 2*n_dim):(5 + 3*n_dim - 1)], # continue integrating qbar
          state[(5 + 3*n_dim):(5 + 4*n_dim - 1)], # continue integrating qbar^2
          state[(5 + 4*n_dim):(5 + 5*n_dim - 1)], # median adaptive now unchanged
          state[(5 + 5*n_dim):(5 + 6*n_dim - 1)] # diagonal elements of S adaptive now unchanged
        )
        
      } else if (whichroot == 2){ # root that leads to adding a new sample to the sample vector used for median estimation
        
        print(t)
        
        transform_qbar_to_q <- state[(5 + 4*n_dim):(5 + 5*n_dim - 1)] + state[(5 + 5*n_dim):(5 + 6*n_dim - 1)] * state[5:(5 + n_dim - 1)] # transform to original q as we are interested in the median of the original q
        
        for(i in 1:n_dim){
          
          update_data_string <- paste("q", i, "_median$add_value(", transform_qbar_to_q[i], ")", sep = "") # add the current values of original q to update estimated median vector
          eval(parse(text = update_data_string))
          
        }
        
        
        newState <- c(
          state[1:3],
          state[4] + time_step_update_data_median, # everything unchanged except for the time step to update the sample vector used to estimate the median vector
          state[5:(5 + 6*n_dim - 1)]
        )
        
        # print(state[4])
        
      } else if (whichroot == 3){
        
        old_m_adaptive <<- m_adaptive
        
        for(i in 1:n_dim){ 
          # at the start of the adaptive period, we update the median right away based on the samples we have collected so far. This is in case we get an event where one of the coordinates crosses the median that was set at initial time. 
          # the elements of S will still be as before since we have not collected any crossing times yet
          update_median_string <- paste("q", i, "_median$q[3]", sep = "") # string that corresponds to the current estimated median of the given coordinate
          
          new_updated_m_adaptive <- eval(parse(text = update_median_string))
          
          m_adaptive[i] <<- new_updated_m_adaptive # set this current estimated median to the median adaptive vector
          
        }
        
        newState <- c(
          state[1:4], # no need to update any of these as these are related to either momentum refresh or updating the sample vector used to estimate the median
          state[5:(5 + n_dim - 1)] + (1 / state[(5 + 5*n_dim):(5 + 6*n_dim - 1)]) * (old_m_adaptive - m_adaptive), # only update qbar now that m has been updated as well. Use this formula since s is unchanged here.
          state[(5 + n_dim):(5 + 4*n_dim - 1)], # pbar and the integrated quantities of qbar are unchanged
          m_adaptive, # median adaptive now has changed
          state[(5 + 5*n_dim):(5 + 6*n_dim - 1)] # diagonal elements of S adaptive unchanged
        )
        
      } else { 
        # register some quantities when the coordinate has crossed the corresponding median, i.e. if qbar has crossed 0
        # also update median adaptive and diagonal elements of S adaptive when this happens
        
        if(t > proportion_time_until_adaptive_start * T){
          
          proposed_s_adaptive <- s_adaptive
          
          proposed_m_adaptive <- m_adaptive
          
          index_q_root <- whichroot - 3
          
          if(n_q_root_median[index_q_root] >= 1){ # only run this part whenever the coordinate has crossed the estimated median at least once as we need two crossing times in order to find a difference between two consecutive times
            
            # New version: Update only the coordinate of m and S related to the coordinate that crosses its median
            # Previously: Only the coordinate of S related to the coordinate that crosses its median was updated. On the other hand, the whole median vector was updated in such events. 
            
            # diff_time_q_root_median[[index_q_root]] <<-  append(diff_time_q_root_median[[index_q_root]], t - cross_median_time[index_q_root]) # store the difference between two consecutive median crossing times
            
            # Dual averaging: Want to find values of diagonal elements of S such that the time between two consecutive events of crossing median is equal to the desired value T*.
            # So if t is the current time when median has been crossed and T_{index_q_root}^{k-1} is the previous time this coordinate crosses median, let alpha_k = t - T_{index_q_root}^{k - 1}.
            # This is equivalent to say that we want to find values of s_{index_q_root} such that the difference between T* and T_{index_q_root}^{k-1} is zero when k gets larger. 
            # Therefore if this difference between T* and T_{index_q_root}^{k-1} is denoted as alpha_k, we can define H_k = T* - alpha_k and use dual averaging. 
            
            update_s_string <- paste("s", index_q_root, "_solution$update_solution(", ((desired_time_until_crossing_median) - (t - cross_median_time[index_q_root])), ")", sep = "") # Add the value of H_k for the current x_k (i.e. current adaptive value of s)
            
            eval(parse(text = update_s_string)) # run the update
            
            update_s_adaptive_string <- paste("s", index_q_root, "_solution$x_bar", sep = "") # the string to get out the new value of s adaptive
            
            new_updated_s_adaptive <- eval(parse(text = update_s_adaptive_string)) # run the evaluation to get out the new value of s adaptive
            
            proposed_s_adaptive[index_q_root] <- exp(new_updated_s_adaptive) # update s adaptive for the coordinate that crosses median at this time point
            
            update_median_string <- paste("q", index_q_root, "_median$q[3]", sep = "") # string that corresponds to the current estimated median of the given coordinate
            
            new_updated_m_adaptive <- eval(parse(text = update_median_string))
            
            proposed_m_adaptive[index_q_root] <- new_updated_m_adaptive # set this current estimated median to the median adaptive vector
            
            old_m_adaptive <<- m_adaptive
            
            old_s_adaptive <<- s_adaptive
            
            # for(i in 1:n_dim){
            #   
            #   update_median_string <- paste("q", i, "_median$q[3]", sep = "") # string that corresponds to the current estimated median of the given coordinate
            #   
            #   new_updated_m_adaptive <- eval(parse(text = update_median_string))
            #   
            #   m_adaptive[i] <<- new_updated_m_adaptive # set this current estimated median to the median adaptive vector
            #   
            # }
            
            m_adaptive <<- proposed_m_adaptive
            
            s_adaptive <<- proposed_s_adaptive
            
            new_qbar <- (1 / s_adaptive) * old_s_adaptive * state[5:(5 + n_dim - 1)] + (1 / s_adaptive) * (old_m_adaptive - m_adaptive) # formula to change qbar when m and S change in order to keep original q unchanged at each update of m and S
            
            newState <- c(
              
              state[1:4], # no need to update any of these as these are related to either momentum refresh or updating the sample vector used to estimate the median
              new_qbar, # qbar needs to be changed in order for q to stay still after an event where m and S have changed
              state[(5 + n_dim):(5 + 2*n_dim - 1)], # pbar unchanged
              state[(5 + 2*n_dim):(5 + 3*n_dim - 1)], # continue integrating qbar
              state[(5 + 3*n_dim):(5 + 4*n_dim - 1)], # continue integrating qbar^2
              m_adaptive, # update median vector
              s_adaptive # update the diagonal elements of S
            )
            
          } else {
            
            newState <- c(
              state[1:(5 + 6*n_dim - 1)]
            )
            
          }
          
          cross_median_time[index_q_root] <<- t # register that current time point is the most recent time point that this specific coordinate has crossed median
          
          # time_q_root_median[[index_q_root]] <<- append(time_q_root_median[[index_q_root]], t) # store the current time point that this specific coordinate has crossed median
          
          n_q_root_median[index_q_root] <<- n_q_root_median[index_q_root] + 1 # update the count of how many times the coordinate has crossed median
          
          
        } else { # keep everything unchanged when crossing median after the adaptive period is over
          
          newState <- c(
            state[1:(5 + 6*n_dim - 1)]
          )
          
        }
        
      }
      
    } else { # in case we get a root at the start, everything is unchanged
      
      newState <- c(
        state[1:(5 + 6*n_dim - 1)]
      )
      
    }
    
    return(newState)
    
  }
  
  if(is.null(pbar_initial)){
    pbar_initial <- rnorm(n_dim)  
  } else {
    pbar_initial <- pbar_initial
  }
  
  u_initial <- rexp(1)
  
  y0 <- c(
    number_of_events = 0,
    Lambda = 0,
    u = u_initial,
    time_update_data_median = time_step_update_data_median,
    qbar = qbar_initial,
    pbar = pbar_initial,
    int_qbar = rep(0, n_dim),
    int_qbar_squared = rep(0, n_dim),
    m = m_adaptive,
    s = s_adaptive
  )
  
  # run actual simulation
  
  sim.out <- deSolve::lsodar(
    y = y0,
    func = ode,
    times = seq(from = 0, to = T, length.out = n_samples + 1), # the first period of warmup corresponds to N*samples + 1
    rootfunc = root.fun, 
    events = list(func=at.event, root = TRUE),
    rtol = rtol,
    atol = atol,
    maxsteps = maxsteps # maximal number of steps per output interval taken by the solver (might need to increase this in certain cases)
  )
  
  # Convert the results into a data frame
  df.sim.out <- as.data.frame(sim.out)
  
  return(
    list(
      output_from_ode_solver = df.sim.out, 
      n_evals_ode = n_evals_ode, 
      m_adaptive = as.numeric(m_adaptive), 
      s_adaptive = as.numeric(s_adaptive),
      qbar_end = as.numeric(df.sim.out[nrow(df.sim.out), 5:(5 + n_dim - 1) + 1]), # plus one since lsodar adds a time column as well
      pbar_end = as.numeric(df.sim.out[nrow(df.sim.out), (5 + n_dim):(5 + 2*n_dim - 1) + 1]),
      u_end = df.sim.out[nrow(df.sim.out), 3 + 1], # value of u at the time the adaptive process ends
      Lambda_end = df.sim.out[nrow(df.sim.out), 2 + 1], # value of Lambda at the time the adaptive process ends
      n_q_root_median = n_q_root_median,
      cross_median_time = cross_median_time
    )
  )
  
  
}


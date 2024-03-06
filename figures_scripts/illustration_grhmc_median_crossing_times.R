# Function to simulate Hamiltonian dynamics using GRHMC with no transformation - used to produce Figure 1

illustration_grhmc_median_crossing_times <- function(log_target_grad,
                                                     lambda,
                                                     T = 10000,
                                                     n_samples = 2000,
                                                     s_elements = c(1, 1), 
                                                     median_elements = c(0, 0),
                                                     q_initial = NULL,
                                                     random_state = NULL,
                                                     maxsteps = NULL,
                                                     rtol = NULL,
                                                     atol = NULL,
                                                     time_step_update_data_median = 1
){
  
  if(!is.null(random_state)){
    set.seed(random_state)
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
  
  n_dim <- length(log_target_grad(0))
  
  n_evals_ode <- 0L # store how many times the ode function below is called in the process
  
  count_n_evals_ode <- function(){
    
    n_evals_ode <<- n_evals_ode + 1 # increment the counter every time the ode function is called
    
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
      s_elements * 
        log_target_grad(median_elements + s_elements * state[5:(5 + n_dim - 1)]), # dot pbar = gradient of log transformed density wrt. qbar
      state[5:(5 + n_dim - 1)], # int qbar dt
      state[5:(5 + n_dim - 1)]^2 # int qbar ^2 dt
    )
    
    return(list(ret))
    
  }
  
  root.fun <- function(t, state, parms){
    
    return(
      
      c(
        state[2] - state[3], # root that leads to standard momentum refresh
        t - state[4], # root that leads to adding a new sample to the sample vector used for median estimation
        state[5:(5 + n_dim - 1)] # events whenever a coordinate crosses its corresponding median value --> qbar crosses 0
      )
      
    )
    
  }
  
  n_q_root_median <- rep(0, n_dim) # store number of times a given position coordinate has crossed median value
  
  diff_time_q_root_median <- rep_len(list(numeric()), n_dim) # store the time between two events of crossing median for each position coordinate
  
  time_q_root_median <- rep_len(list(numeric()), n_dim) # store the specific time points where the median of a given position coordinate has been crossed
  
  cross_median_time <- rep(0, n_dim) # store the fictitious time point when the median is crossed, this will be used to find the time between such events
  
  q_original <- rep_len(list(numeric()), n_dim) # store the specific time points where the median of a given position coordinate has been crossed
  
  at.event <- function(t, state, parms){
    
    yroot <- c(
      state[2] - state[3], # root that leads to standard momentum refresh
      t - state[4], # root that leads to adding a new sample to the sample vector used for median estimation
      state[5:(5 + n_dim - 1)] # events whenever a coordinate crosses its corresponding median value --> qbar crosses 0
    )
    
    whichroot <- which(min(abs(yroot)) == abs(yroot))
    
    if(t > 0){ # when q_initial is equal to initialized median, we get a root right at the start. To avoid this, we set this condition here. 
      
      if(whichroot == 1){ # standard momentum refresh event
        
        p_refresh <- rnorm(n_dim)
        u_refresh <- rexp(1)
        
        newState <- c(
          
          state[1] + 1, # update the number of momentum refresh event counter
          0, # reset Lambda
          u_refresh, # simulate a new u from exponential distribution with rate 1
          state[4], # time to update data for median estimation
          state[5:(5 + n_dim - 1)], # qbar needs to be changed in order for q to stay still after an event where m has changed
          p_refresh, # refresh momentum
          state[(5 + 2*n_dim):(5 + 3*n_dim - 1)], # continue integrating qbar
          state[(5 + 3*n_dim):(5 + 4*n_dim - 1)] # continue integrating qbar^2
          
        )
        
      } else if (whichroot == 2){ # root that leads to adding a new sample to the sample vector used for median estimation
        
        transform_qbar_to_q <- median_elements + s_elements * state[5:(5 + n_dim - 1)] # transform to original q as we are interested in the median of the original q
        
        newState <- c(
          state[1:3],
          state[4] + time_step_update_data_median, # everything unchanged except for the time step to update the sample vector used to estimate the median vector
          state[5:(5 + 4*n_dim - 1)]
        )
        
      } else { # register some quantities when the coordinate has crossed the corresponding median, i.e. if qbar has crossed 0
        
        
        index_q_root <- whichroot - 2
        
        if(n_q_root_median[index_q_root] >= 1){
          
          diff_time_q_root_median[[index_q_root]] <<-  append(diff_time_q_root_median[[index_q_root]], t - cross_median_time[index_q_root])
          
        }
        
        time_q_root_median[[index_q_root]] <<- append(time_q_root_median[[index_q_root]], t)
        
        cross_median_time[index_q_root] <<- t
        
        n_q_root_median[index_q_root] <<- n_q_root_median[index_q_root] + 1
        
        newState <- c(
          state[1:(5 + 4*n_dim - 1)]
        )
        
      }
      
    } else {
      
      newState <- c(
        state[1:(5 + 4*n_dim - 1)]
      )
      
    }
    
    return(newState)
    
  }
  
  if(is.null(q_initial)){
    q_initial <- rep(0, n_dim)
  } else {
    q_initial <- q_initial
  }
  
  p_initial <- rnorm(n_dim)
  u_initial <- rexp(1)
  
  y0 <- c(
    number_of_events = 0,
    Lambda = 0,
    u = u_initial,
    time_update_data_median = time_step_update_data_median,
    qbar = q_initial,
    pbar = p_initial,
    int_q = rep(0, n_dim),
    int_q_squared = rep(0, n_dim)
  )
  
  # run actual simulation
  
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
  
  # Collect samples of the original q after the warmup period
  
  q_samples <- t(median_elements + diag(s_elements) %*% t(df.sim.out[, (5 + 1):(5 + n_dim - 1 + 1)])) # note that the column index is shifted by one here because a time column is added automatically by the ode solver. 
  
  return(
    list(
      q_samples = q_samples, 
      output_from_ode_solver = df.sim.out, 
      n_evals_ode = n_evals_ode, 
      diff_time_q_root_median = diff_time_q_root_median,
      time_q_root_median = time_q_root_median, 
      n_q_root_median = n_q_root_median
    )
  )
  
  
}
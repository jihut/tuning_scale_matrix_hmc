# First step - Tuning m and S using VARI

first_step_grhmc_transformed_function_VARI <- function(log_target_grad, 
                                                       lambda, 
                                                       n_parameters = NULL,
                                                       T = 5000, 
                                                       n_samples = 1000,
                                                       m_initial = NULL,
                                                       diag_s_elements_initial = NULL,
                                                       qbar_initial = NULL,
                                                       pbar_initial = NULL,
                                                       random_state = NULL,
                                                       maxsteps = NULL,
                                                       rtol = NULL,
                                                       atol = NULL,
                                                       proportion_time_until_adaptive_start = 0.05,
                                                       min_proportion_of_previous_state = 0.5,
                                                       max_proportion_of_previous_state = 2
                                                       
){
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  n_dim_error <- FALSE
  
  # define the (first order) ode for state
  # state[1]: number of events up to time point
  # state[2]: Lambda, which resets after each event
  # state[3]: u, simulate from exp(1) after each event
  # state[4:(4 + n_dim - 1)]: qbar
  # state[(4 + n_dim):(4 + 2*n_dim - 1)]: pbar
  # state[(4 + 2*n_dim):(4 + 3*n_dim - 1)]: \int qbar dt
  # state[(4 + 3*n_dim):(4 + 4*n_dim - 1)]: \int qbar^2 dt
  # state[(4 + 4*n_dim):(4 + 5*n_dim - 1)]: \int q dt = \int (m + Sqbar) dt
  # state[(4 + 5*n_dim):(4 + 6*n_dim - 1)]: \int q^2 dt = \int (m + Sqbar)^2 dt
  # state[(4 + 6*n_dim):(4 + 7*n_dim - 1)]: m - mean vector where elements should reflect the mean of the original position coordinates
  # state[(4 + 7*n_dim):(4 + 8*n_dim - 1)]: S - diagonal matrix where the diagonal elements should reflect the marginal standard deviation of the original position coordinates
  
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
  
  ode <- function(t, state, parms){
    
    count_n_evals_ode()
    
    ret <- c(
      
      0, # number of events
      lambda, # integrate lambda to get Lambda
      0, # u
      state[(4 + n_dim):(4 + 2*n_dim - 1)], # \dot qbar = pbar
      state[(4 + 7*n_dim):(4 + 8*n_dim - 1)] * 
        log_target_grad(state[(4 + 6*n_dim):(4 + 7*n_dim - 1)] + state[(4 + 7*n_dim):(4 + 8*n_dim - 1)] * state[4:(4 + n_dim - 1)]), #\dot pbar = gradient of log transformed density wrt. qbar
      state[4:(4 + n_dim - 1)], # \int qbar dt
      state[4:(4 + n_dim - 1)]^2, # \int qbar^2 dt
      state[(4 + 6*n_dim):(4 + 7*n_dim - 1)] + state[(4 + 7*n_dim):(4 + 8*n_dim - 1)] * state[4:(4 + n_dim - 1)], # \int q dt = \int m + Sqbar dt
      (state[(4 + 6*n_dim):(4 + 7*n_dim - 1)] + state[(4 + 7*n_dim):(4 + 8*n_dim - 1)] * state[4:(4 + n_dim - 1)])^2, # \int q^2 dt = \int (m + Sqbar) dt
      rep(0, n_dim), # m - only update at event time
      rep(0, n_dim) # diagonal elements of S - only update at event time
    )
    
    return(list(ret))
    
  }
  
  # events occur when this function evaluate to 0: here \Lambda=u
  
  root.fun <- function(t, state, parms){return(state[2] - state[3])}
  
  # what happens at events?
  
  # Define vectors that changes at each event time in warmup period, representing m and S.
  # Use the final estimates of these two in order to transform back to samples of q original. 
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
  
  at.event <- function(t, state, parms){
    
    if(t > proportion_time_until_adaptive_start * T){
      
      old_m_adaptive <<- m_adaptive
      old_s_adaptive <<- s_adaptive
      
      m_adaptive <<- state[(4 + 4*n_dim):(4 + 5*n_dim - 1)]/t
      proposed_s_adaptive <- sqrt(state[(4 + 5*n_dim):(4 + 6*n_dim - 1)]/t - m_adaptive^2)
      
      s_adaptive <<- pmin(max_proportion_of_previous_state * s_adaptive, pmax(min_proportion_of_previous_state * s_adaptive, proposed_s_adaptive))
      
      pbar_refresh <- rnorm(n_dim)
      u_refresh <- rexp(1)
      
      newState <- c(
        
        state[1] + 1, # update the number of event counter
        0, # reset Lambda
        u_refresh, # simulate a new u from exponential distribution with rate 1
        (1 / s_adaptive) * old_s_adaptive * state[4:(4 + n_dim - 1)] + (1 / s_adaptive) * (old_m_adaptive - m_adaptive), # qbar needs to be changed in order for q to stay still after an event where m and S have both changed
        # state[4:(4 + n_dim - 1)], # qbar unchanged
        pbar_refresh, # p refresh - phi = 0 here so can simulate n_dim iid samples from N(0, 1)
        state[(4 + 2*n_dim):(4 + 6*n_dim - 1)], # continue integrating all the quantities defined 
        m_adaptive, # update m with the new temporal estimates of m
        s_adaptive # update S with the new temporal estimates of the diagonal elements of S. 
        
      )
      
    } else {
      
      pbar_refresh <- rnorm(n_dim)
      u_refresh <- rexp(1)
      
      newState <- c(
        
        state[1] + 1, # update the number of event counter
        0, # reset Lambda
        u_refresh, # simulate a new u from exponential distribution with rate 1
        state[4:(4 + n_dim - 1)], # qbar unchanged
        pbar_refresh, # p refresh - phi = 0 here so can simulate n_dim iid samples from N(0, 1)
        state[(4 + 2*n_dim):(4 + 6*n_dim - 1)], # continue integrating all the quantities defined 
        m_adaptive, # update m with the new temporal estimates of m
        s_adaptive # update S with the new temporal estimates of the diagonal elements of S. 
        
      )
      
    }
    
    return(newState)
    
  }
  
  # initial state
  
  if(is.null(qbar_initial)){
    # qbar_initial = rep(0, n_dim)
    qbar_initial <- rnorm(n_dim) # If initial values of qbar are not specified --> simulate randomly from a multivariate normal distribution with zero mean vector and diagonal covariance matrix
  } else {
    qbar_initial = qbar_initial
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
    qbar = qbar_initial,
    pbar = pbar_initial,
    int_qbar = rep(0, n_dim),
    int_qbar_squared = rep(0, n_dim),
    int_q = rep(0, n_dim),
    int_q_squared = rep(0, n_dim),
    m = m_adaptive,
    s = s_adaptive
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
  
  df.sim.out <- as.data.frame(sim.out)
  
  
  return(
    list(
      output_from_ode_solver = df.sim.out, 
      n_evals_ode = n_evals_ode, 
      m_adaptive = as.numeric(m_adaptive), 
      s_adaptive = as.numeric(s_adaptive),
      qbar_end = as.numeric(df.sim.out[nrow(df.sim.out), 4:(4 + n_dim - 1) + 1]), # plus one since lsodar adds a time column as well
      pbar_end = as.numeric(df.sim.out[nrow(df.sim.out), (4 + n_dim):(4 + 2*n_dim - 1) + 1]),
      u_end = df.sim.out[nrow(df.sim.out), 3 + 1], # value of u at the time the adaptive process ends
      Lambda_end = df.sim.out[nrow(df.sim.out), 2 + 1] # value of Lambda at the time the adaptive process ends
    )
  )
  
}

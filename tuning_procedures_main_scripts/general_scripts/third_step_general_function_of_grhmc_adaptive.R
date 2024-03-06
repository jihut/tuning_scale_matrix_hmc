# Third step in the comparison process - fixed m, S and lambda found from previous two steps and generate samples using GRHMC framework

third_step_grhmc_transformed_function <- function(
    log_target_grad,
    lambda,
    T = 5000,
    n_samples = 1000,
    diag_s_elements_initial,
    m_initial,
    qbar_initial,
    pbar_initial,
    Lambda_initial,
    u_initial,
    random_state = NULL,
    maxsteps = NULL,
    rtol = NULL,
    atol = NULL
){
  
  if(!is.null(random_state)){
    set.seed(random_state)
  }
  
  n_dim <- length(qbar_initial)
  
  n_evals_ode <- 0L # store how many times the ode function below is called in the process
  
  count_n_evals_ode <- function(){
    
    n_evals_ode <<- n_evals_ode + 1 # increment the counter every time the ode function is called
    
  }
  
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
  
  ode <- function(t, state, parms){
    
    count_n_evals_ode()
    
    ret <- c(
      
      0, # number of events
      lambda, # integrate lambda to get Lambda
      0, # u
      state[(4 + n_dim):(4 + 2*n_dim - 1)], # \dot qbar = pbar
      diag_s_elements_initial * 
        log_target_grad(m_initial + diag_s_elements_initial * state[4:(4 + n_dim - 1)]), #\dot pbar = gradient of log transformed density wrt. qbar
      state[4:(4 + n_dim - 1)], # \int qbar dt
      state[4:(4 + n_dim - 1)]^2, # \int qbar^2 dt
      m_initial + diag_s_elements_initial* state[4:(4 + n_dim - 1)], # \int q dt = \int m + Sqbar dt
      (m_initial + diag_s_elements_initial * state[4:(4 + n_dim - 1)])^2 # \int q^2 dt = \int (m + Sqbar) dt
      
    )
    
    return(list(ret))
    
  }
  
  # events occur when this function evaluate to 0: here \Lambda=u
  
  root.fun <- function(t, state, parms){return(state[2] - state[3])}
  
  at.event <- function(t, state, parms){
    
    pbar_refresh <- rnorm(n_dim)
    u_refresh <- rexp(1)
    
    newState <- c(
      
      state[1] + 1, # update the number of event counter
      0, # reset Lambda
      u_refresh, # simulate a new u from exponential distribution with rate 1
      state[4:(4 + n_dim - 1)], # qbar unchanged
      pbar_refresh, # p refresh - phi = 0 here so can simulate n_dim iid samples from N(0, 1)
      state[(4 + 2*n_dim):(4 + 6*n_dim - 1)] # continue integrating all the quantities defined
      
    )
    
    return(newState)
    
  }
  
  y0 <- c(
    number_of_events = 0,
    Lambda = Lambda_initial,
    u = u_initial,
    qbar = qbar_initial,
    pbar = pbar_initial,
    int_qbar = rep(0, n_dim),
    int_qbar_squared = rep(0, n_dim),
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
  
  q_original_samples <- t(m_initial + diag(diag_s_elements_initial) %*% t(df.sim.out[, 4:(4 + n_dim - 1) + 1]))[-1, ]
  
  return(
    list(
      q_original_samples = q_original_samples,
      output_from_ode_solver = df.sim.out,
      n_evals_ode = n_evals_ode,
      s_elements = diag_s_elements_initial,
      m_elements = m_initial,
      lambda = lambda,
      random_state = random_state,
      maxsteps = maxsteps,
      rtol = rtol, 
      atol = atol
    )
  )
  
}

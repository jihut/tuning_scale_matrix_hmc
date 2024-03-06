# Function to tune lambda, the intensity in which momentum refresh occurs

second_step_lambda_adapt_function <- function(fun,q0,
                                              p0,
                                              m.c=rep(0.0,length(q0)),S.c=rep(1.0,length(q0)),
                                              Tmax=5000.0,lambda0=0.2,ema.beta=0.99,fac=1.0,
                                              nsample=10000L,
                                              max.nut.time=100.0,
                                              maxsteps = NULL,
                                              rtol = NULL,
                                              atol = NULL
){
  
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
  
  d <- length(q0)
  lambda <- lambda0
  lambdas <- c(lambda0)
  
  current.refresh.event.state <- q0
  next.refresh.event.state <- q0
  refresh.states <- q0 ##
  nut.found <- FALSE
  refresh.times <- c(0.0)
  nut.times <- c()
  nut.off.val  <- 1.0
  nuts.m <- log(1.0/lambda0)
  
  last.nut.search <- max.nut.time
  
  
  ode.fun <- function(t,z,parms){
    qb <- z[1:d]
    pb <- z[(d+1):(2*d)]
    grad <- fun(m.c + qb*S.c)*S.c
    # return(c(pb,grad,lambda,0.0))
    return(list(c(pb,grad,lambda,0.0)))
  }
  
  
  
  rootFun <- function(t,z,parms){
    #print(nut.found)
    if(nut.found){
      return(c(z[2*d+1]-z[2*d+2],nut.off.val,nut.off.val))
    } else {
      # nut <- sum((z[1:d]-current.refresh.event.state)*z[(d+1):(2:d)])/
      #   sqrt(sum((z[1:d]-current.refresh.event.state)^2)*sum(z[(d+1):(2:d)]^2)+1.0e-8) # this gives error due to 2:d part? fix below by changing to 2*d
      nut <- sum((z[1:d]-current.refresh.event.state)*z[(d+1):(2*d)])/
        sqrt(sum((z[1:d]-current.refresh.event.state)^2)*sum(z[(d+1):(2*d)]^2)+1.0e-8)
      #print(nut)
      return(c(z[2*d+1]-z[2*d+2], # regular momentum refresh
               nut, # NUT-criterion
               t-last.nut.search # stop nut search
      ))
    }
  }
  
  
  
  eventFun <- function(t,z.in,parms){
    if(t>1.0e-6){ # avoid doing anything at the first test run of this function
      # determine which event occurred
      rootFun.at.event <- rootFun(t,z.in,parms)
      event.id <- which.min(abs(rootFun.at.event))
      #print(c(event.id,t,nut.found))
      #print(rootFun.at.event)
      if(event.id==1){
        # momentum refresh
        next.refresh.event.state <<- z.in[1:d]
        refresh.states <<- rbind(refresh.states,z.in[1:d]) ##
        refresh.times <<- c(refresh.times,t)
        
        if(nut.found){
          #print("B1")
          # if nut has been found
          nut.found <<- FALSE
          current.refresh.event.state <<- next.refresh.event.state
          lambda <<- 1.0/(fac*exp(nuts.m))
          lambdas <<- c(lambdas,lambda)
          last.nut.search <<- t+max.nut.time
          return(c(current.refresh.event.state,rnorm(d),0.0,rexp(1)))
        } else {
          #print("B2")
          # otherwise continue integrating
          
          return(c(z.in[1:(2*d+1)],-1.0)) # last effectively turns of momentum refresh search
        }
        
      } else if(event.id==2) {
        # determine in refresh has happened
        if(z.in[2*d+2]>0.0){
          #print("C1")
          # refresh has not yet happened
          nut.found <<- TRUE
          nut.time <- t-refresh.times[length(refresh.times)]
          nuts.m <<- (1.0-ema.beta)*log(nut.time) + ema.beta*nuts.m
          nut.times <<- c(nut.times,nut.time)
          #nut.off.val <<- rootFun.at.event[2]
          return(z.in)
        } else {
          #print("C2")
          # refresh has already happened, step back there and restart
          nut.found <<- FALSE
          nut.time <- t-refresh.times[length(refresh.times)-1]
          nuts.m <<- (1.0-ema.beta)*log(nut.time) + ema.beta*nuts.m
          nut.times <<- c(nut.times,nut.time)
          current.refresh.event.state <<- next.refresh.event.state
          lambda <<- 1.0/(fac*exp(nuts.m))
          lambdas <<- c(lambdas,lambda)
          last.nut.search <<- t+max.nut.time
          return(c(current.refresh.event.state,rnorm(d),0.0,rexp(1)))
        }
      } else {
        print("giving up search for nut event")
        last.nut.search <<- t+max.nut.time
        current.refresh.event.state <<- next.refresh.event.state
        return(c(current.refresh.event.state,rnorm(d),0.0,rexp(1)))
      }
    } else {
      return(z.in)
    }
  }
  
  # z0 <- c(q0,rnorm(d),0.0,rexp(1))
  z0 <- c(q0, p0, 0.0, rexp(1))
  times <- seq(from=0.0,to=Tmax,length.out=nsample+1L)
  out <- deSolve::lsodar(z0,times=times,func=ode.fun,
                         rootfunc=rootFun,
                         events=list(func=eventFun,root=TRUE),
                         rtol = rtol, atol = atol, maxsteps = maxsteps)
  # out <- deSolverRoot(y=z0,times=times,func=ode.fun,rootfunc = rootFun, eventfunc = eventFun)
  
  
  # q.last <- out$samples[nrow(out$samples),2:(d+1)]
  q.last <- out[nrow(out),2:(d+1)]
  p.last <- out[nrow(out), (d+2):(2 * d + 1)]
  
  return(list(ode.out=out,
              # ode.out.events=out$event.samples,
              nut.times=nut.times,
              refresh.times=refresh.times,
              refresh.states=refresh.states, ##
              lambda.last=lambda,
              q.last=q.last,
              p.last = p.last,
              lambdas=lambdas ##
  ))
  
}


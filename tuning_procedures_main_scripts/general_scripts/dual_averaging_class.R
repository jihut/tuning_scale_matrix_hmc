# Script to initialize a class corresponding to the modified dual averaging method by Hoffman et al. (2014) for each position coordinate in the Hamlitonian system
# Using R6 class is faster than reference class. 
library(R6)

dual_averaging_class <- R6Class(
  "dual_averaging_class",
  portable = FALSE,
  cloneable = FALSE,
  class = FALSE,
  public = list(
    
    kappa = NULL,
    gamma = NULL,
    t0 = NULL,
    mu = NULL,
    x = NULL,
    x_bar = NULL,
    sum_H = 0,
    t = 1,
    
    initialize = function(kappa = 0.75, gamma = 0.5, t0 = 10, mu = 1, x_initial){
      kappa <<- kappa
      gamma <<- gamma
      t0 <<- t0
      mu <<- mu
      x <<- x_initial
      x_bar <<- x_initial
    },
    
    eta_t = function(t){
      
      t^(-kappa)
      
    }, 
    
    update_solution = function(H_t){
      
      sum_H <<- sum_H + H_t
      
      eta <- eta_t(t)
      
      x <<- mu - (sqrt(t) / gamma) * (1 / (t + t0)) * sum_H 
      
      x_bar <<- eta * x + (1 - eta) * x_bar
      
      t <<- t + 1
      
    }
    
    
  )
  
)

# Example on usage of class

# linear_func <- function(x){
#   
#   3 * x - 5 + rnorm(1)
#   
# }
# 
# set.seed(42)
# 
# x_initial <- 1
# H_initial <- linear_func(x_initial)
# 
# dual_averaging <- dual_averaging_class$new(x_initial = x_initial)
# dual_averaging$update_solution(H_initial)
# dual_averaging$x
# dual_averaging$x_bar
# dual_averaging$t
# 
# while(dual_averaging$t < 100000){
#   
#   dual_averaging$update_solution(linear_func(dual_averaging$x))
#   
# }
# 
# dual_averaging$t
# dual_averaging$x_bar

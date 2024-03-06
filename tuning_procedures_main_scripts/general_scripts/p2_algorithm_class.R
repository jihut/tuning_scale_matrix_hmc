# Script to initialize a class corresponding to the P^2-algorithm by Jain and Chlamtac (1985) for each position coordinate in the Hamlitonian system
# Using R6 class is faster than reference class. 

library(R6)

p2_algorithm_class <- R6Class(
  "P2_quantile_estimator_algorithm",
  portable = FALSE,
  cloneable = FALSE,
  class = FALSE,
  public = list(
    q = numeric(5),
    ns = numeric(5),
    n = numeric(5),
    count = 0,
    p = 0.5,
    
    
    add_value = function(x) {
      if (count < 5) {
        q[count + 1] <<- x
        
        if (count == 4) {
          q <<- sort(q)
          
          for (i in 1:5) {
            n[i] <<- i
          }
          
          ns[1] <<- 1
          ns[2] <<- 2 * p + 1
          ns[3] <<- 4 * p + 1
          ns[4] <<- 3 + 2 * p
          ns[5] <<- 5
        }
        
        count <<- count + 1
        
        # dns1 <<- 0
        # dns2 <<- p / 2
        # dns3 <<- p
        # dns4 <<- (1 + p) / 2
        # dns5 <<- 1
        
        return()
        
      }
      
      k <- -1 # init
      
      if (x < q[1]) {
        q[1] <<- x
        k <- 1
        
      } else if (x < q[2]) {
        k <- 1
        
      } else if (x < q[3]) {
        k <- 2
        
      } else if (x < q[4]) {
        k <- 3
        
      } else if (x < q[5]) {
        k <- 4
      } else{
        q[5] <<- x
        k <- 4
      }
      
      for (i in (k + 1):5) {
        n[i] <<- n[i] + 1
      }
      
      ns[2] <<- (count * p / 2) + 1
      ns[3] <<- (count * p) + 1
      ns[4] <<-
        (count * (1 + p) / 2) + 1
      ns[5] <<- count + 1
      
      for (i in 2:4) {
        d <- ns[i] - n[i]
        
        if ((d >= 1 &
             ((n[i + 1] - n[i]) > 1)) | (d <= -1 & ((n[i - 1] - n[i]) < -1))) {
          d_int <- as.integer(sign(d))
          qs <- .parabolic(i, d_int)
          
          if (q[i - 1] < qs &
              qs < q[i + 1]) {
            q[i] <<- qs
          } else {
            q[i] <<- .linear(i, d_int)
          }
          n[i] <<- n[i] + d_int
        }
        
        
      }
      
      count <<- count + 1
      
    },
    
    .parabolic = function(i, d) {
      quantity <- q[i] +
        d / (n[i + 1] - n[i - 1]) *
        ((n[i] - n[i - 1] + d) * (q[i + 1] - q[i]) / (n[i +
                                                          1] - n[i]) + (n[i + 1] - n[i] - d) * (q[i] - q[i - 1]) / (n[i] - n[i - 1]))
      
      return(quantity)
      
    },
    
    .linear = function(i, d) {
      quantity <- q[i] + d * (q[i + d] - q[i]) / (n[i + d] - n[i])
      
      return(quantity)
    }
    
  )
  
)


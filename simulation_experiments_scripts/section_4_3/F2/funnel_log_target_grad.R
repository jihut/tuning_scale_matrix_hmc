# Consider a smile shaped distribution in the following form: 

# q1 ~ N(0, 1)
# q2|q1 ~ N(0, exp(2 * q1))

target.grad <- function(q){
  
  return(
    c(
      -1 - q[1] + 1 * exp(-2 * q[1]) * q[2]^2,
      -q[2] * exp(-2 * q[1])
    )
  )
  
}

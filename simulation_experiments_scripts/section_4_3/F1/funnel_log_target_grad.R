# Bivariate funnel distribution - F1 from Table 3

# Consider a smile shaped distribution in the following form: 

# q1 ~ N(0, 1)
# q2|q1 ~ N(0, exp(3 / 2 * q1))

target.grad <- function(q){
  
  return(
    c(
      -3 / 4 - q[1] + 3 / 4 * exp(-3 / 2 * q[1]) * q[2]^2,
      -q[2] * exp(-3 / 2 * q[1])
    )
  )
  
}

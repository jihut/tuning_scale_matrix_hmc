# Bivariate smiley distribution - NG2 from Table 2

# Consider a smile shaped distribution in the following form: 

# q1 ~ N(0, 1)
# q2|q1 ~ N(q1^2, 1)

target.grad <- function(q){
  
  return(
    c(
      -q[1] + 2*q[1] * (q[2] - q[1]^2),
      -(q[2] - q[1]^2)
    )
  )
  
}

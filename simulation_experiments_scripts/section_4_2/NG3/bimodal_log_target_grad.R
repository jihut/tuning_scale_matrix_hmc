# Bivariate bimodal distribution - NG3 from Table 2

# Consider a bivariate bimodal distribution in the following form: 

# Log target density = -(1-q_1^2)^2 – (q_2-q_1)^2/2

target.grad <- function(q){
  
  return(
    c(
      4 * (1 - q[1]^2) * q[1] + (q[2] - q[1]),
      -(q[2] - q[1])
    )
  )
  
}

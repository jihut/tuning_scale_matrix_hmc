# Bivariate t distribution example - NG1 from Table 2

# consider a t distribution with the following setup: 
Sigma <- matrix(
  c(4, 2,
    2, 9),
  nrow = 2,
  byrow = TRUE
)

mu <- c(1, 2)

prec <- solve(Sigma)

# (log) target gradient of the distribution

p <- length(mu) # number of dimensions

v <- 4 # number of degrees of freedoms

target.grad <- function(q) {
  # gradient of log target density based on the given density from wikipedia
  
  return(-c((v + p) / (1 + 1 / v * t(q - mu) %*% prec %*% (q - mu))) * (prec %*% (q - mu)) / v)
  
}

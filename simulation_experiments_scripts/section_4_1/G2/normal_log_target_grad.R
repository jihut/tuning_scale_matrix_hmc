# Normal distribution example - G2 from Table 1

# consider a N(mu, Sigma) target where
Sigma <- matrix(
  c(10, 5,
    5, 10^3),
  nrow = 2,
  byrow = TRUE
)

mu <- c(0, 0)

prec <- solve(Sigma) # precision matrix - inverse of sigma

# target gradient (change to get another target distribution)
target.grad <- function(q) {return(-prec%*%(q - mu))}

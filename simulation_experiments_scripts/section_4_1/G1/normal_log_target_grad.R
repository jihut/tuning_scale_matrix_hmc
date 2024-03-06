# Normal distribution example - G1 from Table 1

# consider a N(mu, Sigma) target where
Sigma <- matrix(
  c(4, 0.5,
    0.5, 9),
  nrow = 2,
  byrow = TRUE
)

mu <- c(1, 2)

prec <- solve(Sigma) # precision matrix - inverse of sigma

# target gradient (change to get another target distribution)
target.grad <- function(q) {return(-prec%*%(q - mu))}

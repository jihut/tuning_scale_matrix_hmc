# Normal distribution example - G4 from Table 1

# consider a N(mu, Sigma) target where
Sigma <- diag(rep(1, 10))

mu <- rep(0, 10)

prec <- solve(Sigma) # precision matrix - inverse of sigma

# target gradient (change to get another target distribution)
target.grad <- function(q) {return(-prec%*%(q - mu))}

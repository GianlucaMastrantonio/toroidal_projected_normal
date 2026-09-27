#cholesky <- function(x, ...) {
#  chol(x, ...)
#}

#dInvWishart <- function(S, nu, Psi, log = TRUE) {
#  LaplacesDemon::dinvwishart(Sigma = S, nu = nu, S = Psi, log = log)
#}

#rInvWishart <- function(n, nu, Psi) {
#  draws <- replicate(n, LaplacesDemon::rinvwishart(nu = nu, S = Psi), simplify = "array")

#  if (n == 1) {
#    dim(draws) <- c(dim(Psi), 1)
#  }

#  draws
#}

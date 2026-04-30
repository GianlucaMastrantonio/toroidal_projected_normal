theta <- 0.3
R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
S <- diag(1, 2)
S[2, 2] <- 1
S[1, 2] <- S[2, 1] <- 1
round(R %*% S %*% t(R), 5)

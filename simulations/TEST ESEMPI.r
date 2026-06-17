theta <- 0.3
R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
S <- diag(1, 2)
S[2, 2] <- 1
S[1, 2] <- S[2, 1] <- 1
round(R %*% S %*% t(R), 5)


par_xi <-
  mu_fc <- 0.5 * (log(2) + digamma(par_xi / 2))
var_fc <- 0.25 * trigamma(par_xi / 2)

dfs <- c(20, 50, 100, 200)

par(mfrow = c(2, 2))

for (nu in dfs) {
  z <- seq(2, 3, length.out = 1000)

  exact <- dchisq(exp(2 * z), df = nu) * 2 * exp(2 * z)

  mu_yours <- 0.5 * log(nu)
  var_yours <- 1 / (2 * nu)

  mu_exact <- 0.5 * (log(2) + digamma(nu / 2))
  var_exact <- 0.25 * trigamma(nu / 2)

  approx_yours <- dnorm(z, mean = mu_yours, sd = sqrt(var_yours))
  approx_mine <- dnorm(z, mean = mu_exact, sd = sqrt(var_exact))

  plot(z, exact,
    type = "l", lwd = 2,
    main = paste("df =", nu),
    xlab = "z = log(sqrt(chi-square))",
    ylab = "density"
  )

  lines(z, approx_yours, lty = 2, lwd = 2)
  lines(z, approx_mine, lty = 3, lwd = 2)

  legend("topright",
    legend = c("exact", "your approx", "moment matched"),
    lty = c(1, 2, 3),
    lwd = 2,
    bty = "n"
  )
}

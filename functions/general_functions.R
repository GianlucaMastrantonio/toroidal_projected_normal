q_wc <- function(u, mu, lambda) {
  # Ensure u is in (0,1)
  if (any(u <= 0 | u >= 1)) {
    stop("u must be strictly between 0 and 1.")
  }

  # Compute z_u
  cos_term <- cos(2 * pi * u)
  numerator <- 2 * lambda + (1 + lambda^2) * cos_term
  denominator <- 1 + lambda^2 + 2 * lambda * cos_term
  z_u <- numerator / denominator

  # Clamp z_u to [-1, 1] to avoid numerical issues
  z_u <- pmin(1, pmax(-1, z_u))

  # Compute quantile
  angle <- acos(z_u)
  a <- ifelse(u <= 0.5, mu + angle, mu - angle)

  # Wrap to [0, 2pi)
  a <- (a %% (2 * pi))

  return(a)
}

cdf_wc_un <- function(theta_un, mu, rho) {
  num <- (1 + rho^2) * cos(theta_un - mu) - 2 * rho
  den <- 1 + rho^2 - 2 * rho * cos(theta_un - mu)
  if (sin(theta_un - mu) >= 0) {
    return((1 / (2 * pi)) * acos(num / den))
  } else {
    return(1 - (1 / (2 * pi)) * acos(num / den))
  }
}
cdf_wc <- function(theta, mu, rho) {
  n <- length(theta)
  ret <- rep(NA, n)
  d0 <- cdf_wc_un(0, mu, rho)
  for (i in 1:n)
  {
    ret[i] <- cdf_wc_un(theta[i], mu, rho)
  }
  # return((ret - d0) %% 1)
  return(ret)
}

func_cdf_wc_un <- function(theta_un, mu, rho) {
  num <- (1 + rho^2) * cos(theta_un - mu) - 2 * rho
  den <- 1 + rho^2 - 2 * rho * cos(theta_un - mu)
  if (sin(theta_un - mu) >= 0) {
    return((1 / (2 * pi)) * acos(num / den))
  } else {
    return(1 - (1 / (2 * pi)) * acos(num / den))
  }
}
func_cdf_wc <- function(theta, mu, rho) {
  n <- length(theta)
  ret <- rep(NA, n)
  d0 <- func_cdf_wc_un(0, mu, rho)
  for (i in 1:n)
  {
    ret[i] <- func_cdf_wc_un(theta[i], mu, rho)
  }
  return((ret - d0) %% 1)
}

func_d_wc <- function(theta, mu, rho) {
  return(1 / (2 * pi) * ((1 - rho^2) / (1 + rho^2 - 2 * rho * cos(theta - mu))))
}

func_logd_wc <- function(theta, mu, rho) {
  return(-log(2 * pi) + log(1 - rho^2) - log(1 + rho^2 - 2 * rho * cos(theta - mu)))
}

sim_sigma <- function(par1, par2) {
  tryCatch(
    {
      Sigma_s <- riwish(par1, par2)
      chol(abs(Sigma_s))
      return(list(Sigma_s, TRUE))
    },
    error = function(e) {
      return(list(1, FALSE))
    }
  )
}
test_sigma <- function(Sigma_s) {
  tryCatch(
    {
      chol(abs(Sigma_s))
      return(TRUE)
    },
    error = function(e) {
      return(FALSE)
    }
  )
}
crps_circ <- function(real_data, missing_vec) {
  dd <- c(real_data, missing_vec)

  dist_mat <- 1 - cos(as.matrix(dist(dd)))
  L <- length(missing_vec)
  return(sum(dist_mat[1, -1]) / L - 1 / (2 * L^2) * sum(c(dist_mat[-1, -1])))
}

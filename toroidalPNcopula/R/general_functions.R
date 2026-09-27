#' Simulate from a Rice distribution
#'
#' Generates Rice-distributed random variables using the Euclidean norm of two
#' independent normal random variables.
#'
#' @param n Number of random values to generate.
#' @param nu Noncentrality parameter.
#' @param sigma Scale parameter. Defaults to `1`.
#'
#' @return A numeric vector of length `n`.
#' @export
#'
#' @examples
#' r_rice(5, nu = 1, sigma = 0.5)
r_rice <- function(n, nu, sigma = 1) {
  x <- rnorm(n, mean = nu, sd = sigma)

  y <- rnorm(n, mean = 0, sd = sigma)

  sqrt(x^2 + y^2)
}

safe_dInvWishart <- function(S, nu, Psi, log = TRUE) {
  out <- tryCatch(
    dInvWishart(S, nu, Psi, log = log),
    error = function(e) if (log) -Inf else 0
  )
  if (length(out) != 1 || is.na(out) || !is.finite(out)) {
    return(if (log) -Inf else 0)
  }
  out
}
logsumexp2 <- function(a, b) {
  if (is.infinite(a) && a == -Inf) {
    return(b)
  }

  if (is.infinite(b) && b == -Inf) {
    return(a)
  }

  m <- max(a, b)

  m + log(exp(a - m) + exp(b - m))
}
#' Wrapped Cauchy quantile function
#'
#' Computes the quantile function of the wrapped Cauchy distribution.
#'
#' @param u Numeric vector of probabilities. Values must be strictly between
#'   `0` and `1`.
#' @param mu Location parameter, in radians.
#' @param lambda Wrapped Cauchy concentration parameter.
#'
#' @return A numeric vector of angles wrapped to `[0, 2 * pi)`.
#' @export
#'
#' @examples
#' q_wc(c(0.25, 0.5, 0.75), mu = pi, lambda = 0.5)
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
# func_cdf_wc <- function(theta, mu, rho) {
#  n <- length(theta)
#  ret <- rep(NA, n)
#  d0 <- func_cdf_wc_un(0, mu, rho)
#  for (i in 1:n)
#  {
#    ret[i] <- func_cdf_wc_un(theta[i], mu, rho)
#  }
#  return((ret - d0) %% 1)
# }
#' Wrapped Cauchy distribution function
#'
#' Computes the distribution function of the wrapped Cauchy distribution.
#'
#' @param theta Numeric vector of angles, in radians.
#' @param mu Location parameter, in radians.
#' @param rho Wrapped Cauchy parameter.
#'
#' @return A numeric vector of cumulative probabilities.
#' @export
#'
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 5)
#' func_cdf_wc(theta, mu = pi, rho = 0.5)
func_cdf_wc <- function(theta, mu, rho) {
  ang <- theta - mu
  z <- ((1 + rho^2) * cos(ang) - 2 * rho) / (1 + rho^2 - 2 * rho * cos(ang))
  z <- pmin(1, pmax(-1, z))

  val <- acos(z) / (2 * pi)
  ret <- ifelse(sin(ang) >= 0, val, 1 - val)

  ang0 <- -mu
  z0 <- ((1 + rho^2) * cos(ang0) - 2 * rho) / (1 + rho^2 - 2 * rho * cos(ang0))
  z0 <- pmin(1, pmax(-1, z0))

  val0 <- acos(z0) / (2 * pi)
  d0 <- ifelse(sin(ang0) >= 0, val0, 1 - val0)

  (ret - d0) %% 1
}

#' Wrapped Cauchy density
#'
#' Computes the wrapped Cauchy density.
#'
#' @param theta Numeric vector of angles, in radians.
#' @param mu Location parameter, in radians.
#' @param rho Wrapped Cauchy parameter.
#'
#' @return A numeric vector of density values.
#' @export
#'
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 5)
#' func_d_wc(theta, mu = pi, rho = 0.5)
func_d_wc <- function(theta, mu, rho) {
  return(1 / (2 * pi) * ((1 - rho^2) / (1 + rho^2 - 2 * rho * cos(theta - mu))))
}

#' Wrapped Cauchy log-density
#'
#' Computes the logarithm of the wrapped Cauchy density.
#'
#' @param theta Numeric vector of angles, in radians.
#' @param mu Location parameter, in radians.
#' @param rho Wrapped Cauchy parameter.
#'
#' @return A numeric vector of log-density values.
#' @export
#'
#' @examples
#' theta <- seq(0, 2 * pi, length.out = 5)
#' func_logd_wc(theta, mu = pi, rho = 0.5)
func_logd_wc <- function(theta, mu, rho) {
  return(-log(2 * pi) + log(1 - rho^2) - log(1 + rho^2 - 2 * rho * cos(theta - mu)))
}

#' Simulate a valid covariance matrix
#'
#' Draws a covariance matrix from an inverse-Wishart distribution and checks
#' whether the matrix remains positive definite after taking absolute values.
#'
#' @param par1 Degrees-of-freedom parameter passed to the inverse-Wishart
#'   sampler.
#' @param par2 Scale matrix passed to the inverse-Wishart sampler.
#'
#' @return A list with two elements: the simulated matrix, or `1` on failure,
#'   and a logical flag indicating whether the Cholesky check succeeded.
#' @export
sim_sigma <- function(par1, par2) {
  tryCatch(
    {
      Sigma_s <- rInvWishart(1, par1, par2)[, , 1]
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
test_sigma_mcmc <- function(Sigma_s) {
  tryCatch(
    {
      c1 <- chol(abs(Sigma_s))
      c2 <- chol(Sigma_s)

      return(list(chol_sigma_s = c2, chol_sigma_c = c1, ind = TRUE))
    },
    error = function(e) {
      return(list(ind = FALSE))
    }
  )
}
#' Circular CRPS
#'
#' Computes a circular continuous ranked probability score for an observed
#' angular value and posterior draws or imputations.
#'
#' @param real_data Observed angular value.
#' @param missing_vec Numeric vector of posterior draws or imputations.
#'
#' @return A numeric score.
#' @export
#'
#' @examples
#' crps_circ(real_data = 1, missing_vec = c(0.8, 1.1, 1.3))
crps_circ <- function(real_data, missing_vec) {
  dd <- c(real_data, missing_vec)

  dist_mat <- 1 - cos(as.matrix(dist(dd)))
  L <- length(missing_vec)
  return(sum(dist_mat[1, -1]) / L - 1 / (2 * L^2) * sum(c(dist_mat[-1, -1])))
}

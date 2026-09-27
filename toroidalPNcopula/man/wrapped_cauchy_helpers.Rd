\name{wrapped_cauchy_helpers}
\alias{q_wc}
\alias{func_cdf_wc}
\alias{func_d_wc}
\alias{func_logd_wc}
\title{Wrapped Cauchy helper functions}
\usage{
q_wc(u, mu, lambda)
func_cdf_wc(theta, mu, rho)
func_d_wc(theta, mu, rho)
func_logd_wc(theta, mu, rho)
}
\arguments{
\item{u}{Probability value in the open interval \eqn{(0, 1)}.}
\item{theta}{Angular value.}
\item{mu}{Location parameter.}
\item{lambda, rho}{Wrapped Cauchy parameter.}
}
\value{Numeric vector.}
\description{Quantile, distribution, density, and log-density helpers for the wrapped Cauchy distribution.}

\name{misc_helpers}
\alias{r_rice}
\alias{crps_circ}
\title{Miscellaneous circular-data helpers}
\usage{
r_rice(n, nu, sigma = 1)
crps_circ(real_data, missing_vec)
}
\arguments{
\item{n}{Number of observations.}
\item{nu, sigma}{Rice-distribution parameters.}
\item{real_data}{Observed angular value.}
\item{missing_vec}{Posterior draws or imputations for a missing angular value.}
}
\value{Numeric value or vector.}
\description{Utility functions used by the toroidal model samplers and analyses.}

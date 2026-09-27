\name{mcmc_tpn}
\alias{mcmc_tpn}
\title{Fit a toroidal projected normal model by MCMC}
\usage{
mcmc_tpn(theta, burnin, thin, iterations, prior_mu_mean, prior_mu_var,
  prior_kappa_mean, prior_kappa_var, prior_sigma_nu, prior_sigma_psi,
  mu_init, kappa_init, sigma_init, r_init, adapt_batch, adapt_a,
  adapt_b, adapt_alpha_target, sd_mu_scal, par_sigma_adapt,
  na_index = list(NA), do_only_ESS = TRUE, n_test_sigma = 10,
  type_ess = 1, do_ind = FALSE)
}
\arguments{
\item{theta}{Matrix of angular observations.}
\item{burnin, thin, iterations}{MCMC control parameters.}
\item{prior_mu_mean, prior_mu_var, prior_kappa_mean, prior_kappa_var, prior_sigma_nu, prior_sigma_psi}{Prior parameters.}
\item{mu_init, kappa_init, sigma_init, r_init}{Initial values.}
\item{adapt_batch, adapt_a, adapt_b, adapt_alpha_target, sd_mu_scal, par_sigma_adapt}{Adaptive sampler controls.}
\item{na_index}{Optional missing-value indices.}
\item{do_only_ESS, n_test_sigma, type_ess, do_ind}{Sampler options.}
}
\value{A list containing posterior samples and diagnostics.}
\description{Runs the MCMC sampler for the toroidal projected normal model.}

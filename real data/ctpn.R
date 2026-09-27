# Real-data run for the wrapped Cauchy copula model.
# This file is called by real data/1 - launch_ctpn.R with an args vector that
# selects the seed, covariance structure, sampler settings, and model variant.

library(CholWishart)

library(ggplot2)
library(dplyr)
library(readr)
library(lubridate)
library(VGAM)
library(MCMCpack)
library(Rfast)
library(MCMCpack)
library(MASS)
library(truncnorm)
library(matrixcalc)
library(LaplacesDemon)
library(toroidalPNcopula)

# Return the empirical mode of a vector. Used for mixture allocation summaries.
findmode <- function(x) {
    TT <- table(as.vector(x))
    return(as.numeric(names(TT)[TT == max(TT)][1]))
  }

# Load the real angular data and station metadata.
load("real data/data/gauge.Rdata")


# Decode the run selected by the launch script.
seed <- as.integer(args[1])
do_best_init = c(T,F)[as.integer(args[2])]
do_small = c("T1","T2", "T3", "T4", "T5")[as.integer(args[3])]

do_only_ESS = c(TRUE,FALSE)[as.integer(args[4])]
type_ess <- as.integer(args[5])
n_test_sigma <- as.integer(args[6])

do_ind <- c(TRUE,FALSE)[as.integer(args[7])]
mmmolt_iter <- c(1,2)[as.integer(args[8])]

mixt <- c(TRUE,FALSE)[as.integer(args[9])]
Kmax <- as.integer(args[10])

# Include the model and sampler settings in the output prefix.
name_sim <- paste(name_sim, "IND", do_ind,"_", do_only_ESS, type_ess, n_test_sigma, "do_best_init=", do_best_init,"do_small=", do_small, "mmmolt_iter=", mmmolt_iter, "mixt=",mixt, "Kmax=", Kmax, sep = "")
# ========
# * SECTION - data
# ========
set.seed(1)

n <- nrow(theta)
d <- ncol(theta)

# Build the missing-value index used by the sampler. Missing angles are
# initialized uniformly on the circle before fitting.
set.seed(123)
na_list <- list()
theta_no_na <- theta

for (id in 1:d)
{
  na_list[[id]] <- which(is.na(theta_no_na[, id]))
  if (length(na_list[[id]]) > 0) {
    theta_no_na[na_list[[id]], id] <- runif(length(na_list[[id]]), 0, 2 * pi)
  }
}
theta_no_na <- theta_all
n <- nrow(theta_no_na)
d <- ncol(theta_no_na)


# ========
# * SECTION - MCMC
# ========

set.seed(seed)
mmm <-  m_mcmc

# Initialize the MCMC chain. The launch scripts use the deterministic
# data-based initialization.
if(do_best_init == TRUE)
{
  mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x), na.rm=T), sum(cos(x), na.rm=T)))
  rho_init <- rep(1, d)
  r_init <- matrix(1, n, d)
  x_init <- matrix(0, n, d)
  y_init <- matrix(0, n, d)
  rho_seq <- seq(0.05, 0.95, by = 0.01)
  for (i in 1:d)
  {
    rrr <- rep(NA, length(rho_seq))
    for(ik in 1:length(rho_seq))
    {
      rrr[ik] <- sum(func_logd_wc(theta[, i], mu_init[i], rho_seq[ik]), na.rm=T)
    }
    rho_init[i] <- rho_seq[which.max(rrr)]
    
    theta_app <- 2 * pi * func_cdf_wc(theta[, i] - mu_init[i], 0, rho_init[i])
    if(sum(is.na(theta_app)) > 0)
    {
      theta_app[is.na(theta_app)] <- runif(sum(is.na(theta_app)), 0, 2 * pi)
    }  
    
    u <- runif(n)
    r_direct <- sqrt(-2 * log(u))
    r_app <- r_direct
    x_app <- r_app * cos(theta_app)
    y_app <- r_app * sin(theta_app)
    

    y_init[, i] <- y_app
  }
  
  if(do_ind == TRUE)
  {
    sigma_init <- diag(1, d)
  }else{
      sigma_init <- cov(y_init)
    while (inherits(try(chol(abs(sigma_init)), silent = TRUE), "try-error")) {
      
      sigma_init <- sigma_init + diag(0.01,d)

    }
  }
  
}else{
  mu_init <- runif(d, 0,2*pi)
  mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x), na.rm=T), sum(cos(x), na.rm=T)))
  rho_init <- runif(d, 0.5,0.9)
  r_init <- matrix(runif(n*d, 0.8,1.2), n, d)
  sigma_init <- diag(1, d)
}
if(mixt == FALSE)
{
    start <- Sys.time()

  # Fit the non-mixture wrapped Cauchy copula model.
  out_mcmc <- mcmc_cwc(
    theta = theta_no_na, # the circualr data
    burnin = burnin_mcmc * mmm*mmmolt_iter, # burnin
    thin = thin_mcmc * mmm*mmmolt_iter, # thin
    iterations = iter_mcmc * mmm*mmmolt_iter, # total interations
    prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
    prior_mu_var = rep(100000, d),
    prior_rho_a = rep(1, d), # the prior for B()
    prior_rho_b = rep(1, d),
    prior_sigma_nu = d + nu_app, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
    prior_sigma_psi = diag(1, d)*(d+nu_app -d -1),



    # this section set the initial values of the parameters
    mu_init = mu_init,
    rho_init = rho_init,
    sigma_init = sigma_init,
    r_init = r_init,

    # parameters for the adaptive part of Metropolis
    adapt_batch = batch_mcmc,
    adapt_a = a_mcmc,
    adapt_b = b_mcmc,
    adapt_alpha_target = alpha_target,
    sd_mu_scal = 1,
    sd_rho_scal = 0.1,
    par_sigma_adapt = par_adapt_mcmc,
    na_index = na_list,
    n_test_sigma = n_test_sigma,
    do_only_ESS = do_only_ESS,
    type_ess = type_ess,
    do_ind = do_ind
  )
  end <- Sys.time()
  runtime <- end - start

  # ========
  # * SECTION - Output
  # ========
  # Extract posterior samples and apply the same identification transform used
  # in the simulation scripts.
  mu_out <- out_mcmc$mu_out
  rho_out <- out_mcmc$rho_out
  sigma_s_out <- out_mcmc$sigma_s_out
  sigma_c_out <- out_mcmc$sigma_c_out
  r_out <- out_mcmc$r_out

  waic <- out_mcmc$waic
  ### identification
  mu_out <- mu_out %% (2 * pi)
  nsim <- nrow(mu_out)

  for (isim in 1:nsim)
  {
    ss <- matrix(sigma_s_out[isim, ], nrow = d)
    B <- diag(1 / diag(ss)^0.5)

    sigma_s_out[isim, ] <- B %*% matrix(sigma_s_out[isim, ], nrow = d) %*% B
    sigma_c_out[isim, ] <- B %*% matrix(sigma_c_out[isim, ], nrow = d) %*% B
  }


  # Save the full workspace so the post-analysis script can access posterior
  # samples, data, settings, and diagnostics.
  save.image(paste("real data/output/",name_sim,"cwc_seed",seed,".Rdata", sep = ""))
}else{
      start <- Sys.time()

  # Optional mixture branch kept for archived experiments.
  out_mcmc <- mcmc_cwc_mixture(
    theta = theta_no_na, # the circualr data
    burnin = burnin_mcmc * mmm*mmmolt_iter, # burnin
    thin = thin_mcmc * mmm*mmmolt_iter, # thin
    iterations = iter_mcmc * mmm*mmmolt_iter, # total interations
    prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
    prior_mu_var = rep(100000, d),
    prior_rho_a = rep(1, d), # the prior for B()
    prior_rho_b = rep(1, d),
    prior_sigma_nu = d + nu_app, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
    prior_sigma_psi = diag(1, d)*(d+nu_app -d -1),



    # this section set the initial values of the parameters
    mu_init = mu_init,
    rho_init = rho_init,
    sigma_init = sigma_init,
    r_init = r_init,

    # parameters for the adaptive part of Metropolis
    adapt_batch = batch_mcmc,
    adapt_a = a_mcmc,
    adapt_b = b_mcmc,
    adapt_alpha_target = alpha_target,
    sd_mu_scal = 1,
    sd_rho_scal = 0.1,
    par_sigma_adapt = par_adapt_mcmc,
    na_index = na_list,
    n_test_sigma = n_test_sigma,
    do_only_ESS = do_only_ESS,
    type_ess = type_ess,
    do_ind = do_ind
  )
  end <- Sys.time()
  runtime <- end - start

  # ========
  # * SECTION - Output
  # ========
  # Extract mixture posterior samples and identify covariance scale within each
  # component.
  mu_out <- out_mcmc$mu_out
  rho_out <- out_mcmc$rho_out
  sigma_s_out <- out_mcmc$sigma_s_out
  sigma_c_out <- out_mcmc$sigma_c_out
  r_out <- out_mcmc$r_out
  z_out <- out_mcmc$z_out
  waic <- out_mcmc$waic
  ### identification
  mu_out <- mu_out %% (2 * pi)
  nsim <- nrow(mu_out)

  for(k in 1:dim(sigma_s_out)[3])
  {
    for (isim in 1:nsim)
    {
      ss <- matrix(sigma_s_out[isim, ,k], nrow = d)
      B <- diag(1 / diag(ss)^0.5)

      sigma_s_out[isim, ,k] <- B %*% matrix(sigma_s_out[isim, ,k], nrow = d) %*% B
      sigma_c_out[isim, ,k] <- B %*% matrix(sigma_c_out[isim, ,k], nrow = d) %*% B
    }
  }
  


  # Save the full workspace so the post-analysis script can access posterior
  # samples, data, settings, and diagnostics.
  save.image(paste("real data/output/",name_sim,"cwc_seed",seed,".Rdata", sep = ""))
}


# Simulation run for one wrapped Cauchy copula scenario.
# This file is called by simulations/1 - launch_ctpn.R with an args vector that
# selects the dimension, data replicate, covariance setting, chain, and sampler
# options for a single run.

library(CholWishart)

library(MCMCpack)
library(MASS)
library(truncnorm)
library(matrixcalc)
library(LaplacesDemon)
library(Rfast)
library(toroidalPNcopula)

#### #### #### #### #### ####
#### Simulation
#### #### #### #### #### ####

# Containers used while looping over selected simulation settings.
out <- list()
parout <- list()
counter <- 1

# Seed lists make the chain seed, data seed, and covariance seed independent.
seed_list_chain <- 1:200
seed_list_data <- 1:999
seed_list_sigma <- 1:999

# Decode the scenario selected by the launch script.
app_d <- as.integer(args[1]) # 1:2
app_rho <- as.integer(args[2]) # 1:4
app_sigma_ind_dep <- as.integer(args[3]) # 1:2
app_chain <- as.integer(args[4]) # 1:20
app_sigma <- as.integer(args[5]) # 1:...
app_data <- as.integer(args[6]) # 1:....
app_n <- as.integer(args[7]) # 1:3
do_best_init <- c(T, F)[as.integer(args[8])]
do_only_ESS <- c(TRUE, FALSE)[as.integer(args[9])]
type_ess <- as.integer(args[10])
n_test_sigma <- as.integer(args[11])

# Include the sampler options in the output prefix.
name_sim <- paste(name_sim, do_only_ESS, type_ess, n_test_sigma, "do_best_init=", do_best_init, sep = "")
print(name_sim)
select_d <- app_d

for (select_n in app_n:app_n)
{
  for (select_rho in app_rho:app_rho)
  {
    for (select_sigma in app_sigma_ind_dep:app_sigma_ind_dep)
    {
      for (select_chain in app_chain:app_chain)
      {
        for (select_data in app_data:app_data)
        {
          for (select_sigma_type in app_sigma:app_sigma)
          {
            seed_chain <- seed_list_chain[app_chain] * 1000
            seed_data <- seed_list_data[select_data]
            seed_sigma <- seed_list_sigma[select_sigma_type]

            # Select the dimension, sample size, and base latent parameters for
            # this simulation scenario.
            d <- c(50, 25, 5, 40)[select_d] # dimension of the torus
            n <- c(d * 7, d * 15, d * 7, d * 15, d * 30, d * 45)[select_n] # number of observation


            dmax <- max(d)
            ### parameters
            mu <- rep(0, d)
            kappa <- rep(0, d)

            # Build the two covariance scenarios. The first is independent, the
            # second has dependence generated from a larger valid covariance
            # matrix and an exponential distance correlation structure.
            ddddd <- 500
            sigma_array <- array(0, c(ddddd, ddddd, 2))
            sigma_array[, , 1] <- diag(1, ddddd)


            set.seed(seed_sigma)
            d_test <- 6
            max_attempts <- 100000
            attempts <- 0
            repeat{
              attempts <- attempts + 1
              Sigma_try <- sim_sigma(d_test + 1, diag(1, d_test))
              if (Sigma_try[[2]] || attempts >= max_attempts) {
                break
              }
            }
            dist_mat <- as.matrix(dist(1:100))
            sigma_dist <- exp(-1.6 * dist_mat)
            Sigma_try[[1]] <- kronecker(sigma_dist, Sigma_try[[1]][1:5, 1:5])

            sigma_array[, , 2] <- Sigma_try[[1]]
            Sigma_s <- sigma_array[1:d, 1:d, select_sigma]
            Sigma_c <- abs(Sigma_s)
            chol(Sigma_c)
            chol(Sigma_s)

            # Rescale covariance matrices to unit marginal variances.
            B <- diag(1 / diag(Sigma_s^0.5))
            Sigma_s <- B %*% Sigma_s %*% B
            Sigma_c <- B %*% Sigma_c %*% B

            set.seed(seed_data)

            # Simulate the latent linear variables used to generate the copula.
            if (select_n <= 2) {
              x_c <- mvrnorm(n, kappa, Sigma_c)
              x_s <- mvrnorm(n, rep(0, d), Sigma_s)
            } else {
              x_c <- mvrnorm(d * 45, kappa, Sigma_c)[1:n, ]
              x_s <- mvrnorm(d * 45, rep(0, d), Sigma_s)[1:n, ]
            }


            # Project the latent variables onto the torus. These angles define
            # the copula scale before applying the wrapped Cauchy marginals.
            theta_cop <- matrix(NA, nrow = n, ncol = d)
            for (iobs in 1:n)
            {
              for (id in 1:d)
              {
                theta_cop[iobs, id] <- (atan2(x_s[iobs, id], x_c[iobs, id]) + mu[id]) %% (2 * pi)
              }
            }
            r <- matrix(NA, nrow = n, ncol = d)
            for (iobs in 1:n)
            {
              for (id in 1:d)
              {
                r[iobs, id] <- (x_c[iobs, id]^2 + x_s[iobs, id]^2)^0.5
              }
            }

            ##### ##### ##### ##### ##### ##### #####
            ##### wrapped cauchy marginals
            ##### ##### ##### ##### ##### ##### #####

            # Select the wrapped Cauchy marginal parameters.
            mu <- rep(c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6), times = 20)[1:d]

            rho_mat <- matrix(NA, ncol = dmax, nrow = 4)
            rho_mat[1, ] <- rep(0.3, dmax) / 1
            rho_mat[2, ] <- rep(0.6, dmax) / 1
            rho_mat[3, ] <- rep(0.9, dmax) / 1
            rho_mat[4, ] <- (rep(c(0.3, 0.6, 0.9), (dmax + 3) / 3) / 1)[1:dmax]
            rho <- rho_mat[select_rho, 1:d]


            # Transform the copula-scale angles to wrapped Cauchy marginals.
            theta_seq <- seq(0, 2 * pi, by = 0.00001)
            cumulative_wrappedcauchy <- list()
            for (id in 1:d)
            {
              cumulative_wrappedcauchy[[id]] <- cdf_wc(theta_seq, 0, rho[id])
            }


            theta <- matrix(NA, nrow = n, ncol = d)
            for (iobs in 1:n)
            {
              for (id in 1:d)
              {
                theta[iobs, id] <- (q_wc(theta_cop[iobs, id] / (2 * pi), 0, rho[id]) + mu[id]) %% (2 * pi)
              }
            }

            #### #### #### #### #### ####
            #### MCMC function
            #### #### #### #### #### ####

            # Store the true parameters and simulated data used in this run.
            par_list <- list(
              theta = theta,
              theta_cop = theta_cop,
              x_c = x_c,
              x_s = x_s,
              r = r,
              n = n,
              d = d,
              mu = mu,
              rho = rho,
              Sigma_s = Sigma_s,
              Sigma_c = Sigma_c
            )

            # Initialize the MCMC chain.
            set.seed(seed_chain)


            if (do_best_init == TRUE) {
              mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x)), sum(cos(x))))
              rho_init <- rep(1, d)
              r_init <- matrix(1, n, d)
              x_init <- matrix(0, n, d)
              y_init <- matrix(0, n, d)
              rho_seq <- seq(0.05, 0.95, by = 0.01)
              for (i in 1:d)
              {
                rrr <- rep(NA, length(rho_seq))
                for (ik in 1:length(rho_seq))
                {
                  rrr[ik] <- sum(func_logd_wc(theta[, i], mu_init[i], rho_seq[ik]))
                }
                rho_init[i] <- rho_seq[which.max(rrr)]

                theta_app <- 2 * pi * func_cdf_wc(theta[, i] - mu_init[i], 0, rho_init[i])
                u <- runif(n)
                r_direct <- sqrt(-2 * log(u))
                r_app <- r_direct
                x_app <- r_app * cos(theta_app)
                y_app <- r_app * sin(theta_app)


                y_init[, i] <- y_app
              }
              sigma_init <- cov(y_init)
              sigma_init <- cov(y_init)
              while (inherits(try(chol(abs(sigma_init)), silent = TRUE), "try-error")) {
                sigma_init <- sigma_init + diag(0.01, d)
              }
            } else {
              mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x)), sum(cos(x))))
              rho_init <- runif(d, 0.5, 0.9)
              r_init <- matrix(runif(n * d, 0.8, 1.2), n, d)
              sigma_init <- diag(1, d)
            }
            mmm <- m_mcmc
            start <- Sys.time()

            # Fit the wrapped Cauchy copula model.
            out_mcmc <- mcmc_cwc(
              theta = theta, # the circualr data
              burnin = burnin_mcmc * mmm, # burnin
              thin = thin_mcmc * mmm, # thin
              iterations = iter_mcmc * mmm, # total interations
              prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
              prior_mu_var = rep(100000, d),
              prior_rho_a = rep(1, d), # the prior for B()
              prior_rho_b = rep(1, d),
              prior_sigma_nu = d + nu_app, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
              prior_sigma_psi = diag(1, d) * (d + nu_app - d - 1),


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
              n_test_sigma = n_test_sigma,
              do_only_ESS = do_only_ESS,
              type_ess = type_ess
            )
            end <- Sys.time()
            runtime <- end - start

            # Extract posterior samples and apply the same identification
            # transformation used for the true parameters.
            mu_out <- out_mcmc$mu_out
            rho_out <- out_mcmc$rho_out
            sigma_s_out <- out_mcmc$sigma_s_out
            sigma_c_out <- out_mcmc$sigma_c_out
            r_out <- out_mcmc$r_out
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

            # Bundle the true values, data, posterior samples, and diagnostics
            # saved by the simulation study.
            res_list <- list(
              "runtime" = runtime,
              "n" = n,
              "d" = d,
              "rho" = rho,
              "mu" = mu,
              "Sigma_s" = Sigma_s,
              "Sigma_c" = Sigma_c,
              "x_c" = x_c,
              "x_s" = x_s,
              "theta" = theta,
              "rho_out" = rho_out,
              "mu_out" = mu_out,
              "sigma_s_out" = sigma_s_out,
              "sigma_c_out" = sigma_c_out,
              "pos_def_sigma" = out_mcmc$ess_acc_sigma
            )

            # Save one result file per scenario and chain.
            save(res_list, file = paste(
              "simulations/output/", name_sim, "cwc_simulations_results -",
              " select_d=", select_d,
              " select_n=", select_n,
              " select_rho=", select_rho,
              " select_sigma=", select_sigma,
              " select_chain=", select_chain,
              " seed_sigma=", select_sigma_type,
              " seed_data=", select_data,
              ".Rdata"
            ))
            counter <- counter + 1
          }
        }
      }
    }
  }
}

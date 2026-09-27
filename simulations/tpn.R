# Simulation run for one Toroidal Projected Normal scenario.
# This file is called by simulations/1 - launch_tpn.R with an args vector that
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
app_d  = as.integer(args[1]) # 1:2
app_k = as.integer(args[2]) # 1:4
app_sigma_ind_dep = as.integer(args[3]) # 1:2
app_chain = as.integer(args[4]) # 1:200
app_sigma = as.integer(args[5]) # 1:...
app_data = as.integer(args[6]) # 1:....
app_n = as.integer(args[7]) # 1:3
do_best_init = c(T,F)[as.integer(args[8])]
do_only_ESS = c(TRUE,FALSE)[as.integer(args[9])]
type_ess <- as.integer(args[10])
n_test_sigma <- as.integer(args[11])

select_d <- app_d

# Include the sampler options in the output prefix.
name_sim <- paste(name_sim,do_only_ESS, type_ess, n_test_sigma, "do_best_init=", do_best_init, sep = "")
print(name_sim)

for (select_n in app_n:app_n)
{
  for (select_kappa in app_k:app_k)
  {
    for (select_sigma in app_sigma_ind_dep:app_sigma_ind_dep)
    {
      for (select_chain in app_chain:app_chain)
      {
        for(select_data in app_data:app_data)
        {
          for(select_sigma_type in app_sigma:app_sigma)
          {
            seed_chain <- seed_list_chain[app_chain]*1000
            seed_data <- seed_list_data[select_data]
            seed_sigma <- seed_list_sigma[select_sigma_type]

            # Select the dimension, sample size, and true model parameters for
            # this simulation scenario.
            d <- c( 50,25,5, 40)[select_d] # dimension of the torus
            n <- c(d*7, d*15, d*7, d*15, d*30, d* 45)[select_n] # number of observation
            
            dmax <- max(d)
            ### parameters
            mu <- rep(c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6), times = 20)[1:d]

            kappa_mat <- matrix(NA, ncol = dmax, nrow = 4)
            kappa_mat[1, ] <- rep(0.49, dmax)
            kappa_mat[2, ] <- rep(1.1, dmax)
            kappa_mat[3, ] <- rep(2.45, dmax)
            kappa_mat[4, ] <- (rep(c(0.49, 1.1, 2.45), (dmax + 3) / 3))[1:dmax]
            kappa <- kappa_mat[select_kappa, 1:d]

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
            dist_mat = as.matrix(dist(1:100))
            sigma_dist <- exp(-1.6*dist_mat)
            Sigma_try[[1]] <- kronecker(sigma_dist, Sigma_try[[1]][1:5,1:5])
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

            # Simulate the latent linear variables.
            if(select_n <= 2)
            {
              x_c <- mvrnorm(n, kappa, Sigma_c)
              x_s <- mvrnorm(n, rep(0, d), Sigma_s)
            }else{
              x_c <- mvrnorm(d* 45, kappa, Sigma_c)[1:n,]
              x_s <- mvrnorm(d* 45, rep(0, d), Sigma_s)[1:n,]
            }
            

            # Project the latent variables onto the torus and compute the
            # latent radii.
            theta <- matrix(NA, nrow = n, ncol = d)
            for (iobs in 1:n)
            {
              for (id in 1:d)
              {
                theta[iobs, id] <- (atan2(x_s[iobs, id], x_c[iobs, id]) + mu[id]) %% (2 * pi)
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

            #### #### #### #### #### ####
            #### MCMC function
            #### #### #### #### #### ####

            # Store the true parameters and simulated data used in this run.
            par_list <- list(
              theta = theta,
              x_c = x_c,
              x_s = x_s,
              r = r,
              n = n,
              d = d,
              mu = mu,
              kappa = kappa,
              Sigma_s = Sigma_s,
              Sigma_c = Sigma_c
            )

            
            # Initialize the MCMC chain.
            set.seed(seed_chain)
            if(do_best_init == TRUE)
            {
              mean_init <- mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x)), sum(cos(x))))
              kappa_init <- rep(1, d)
              r_init <- matrix(1, n, d)
              x_init <- matrix(0, n, d)
              y_init <- matrix(0, n, d)

              for (i in 1:d)
              {
                theta_app <- theta[, i] - mean_init[i]
                u <- runif(n)
                C <- mean(cos(theta))
                S <- mean(sin(theta))
                R <- sqrt(C^2 + S^2)
                V <- 1 - R
                kappa_init[i] <- 1 / sqrt(2 * V)

                r_app <- r_rice(n, kappa_init, sigma = 1)


                x_app <- r_app * cos(theta_app)
                y_app <- r_app * sin(theta_app)
                y_init[, i] <- y_app
              }
              sigma_init <- cov(y_init)
              while (inherits(try(chol(abs(sigma_init)), silent = TRUE), "try-error")) {

                sigma_init <- sigma_init + diag(0.01,d)

              }
            }else{
              mean_init <- mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x)), sum(cos(x))))
              kappa_init <- runif(d,0.7,0.9)
              r_init <- matrix(runif(n*d, 0.8,1.2), n, d)
              sigma_init <- diag(1, d)
            }
            

            mmm <- m_mcmc
            start <- Sys.time()

            # Fit the Toroidal Projected Normal model.
            out_mcmc <- mcmc_tpn(
              theta = theta, # the circualr data
              burnin = burnin_mcmc * mmm, # burnin
              thin = thin_mcmc * mmm, # thin
              iterations = iter_mcmc * mmm, # total interations
              prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
              prior_mu_var = rep(100000, d),
              prior_kappa_mean = matrix(0, nrow = d, ncol = 1), # the prior for k is  TN(prior_kappa_mean,prior_kappa_var )
              prior_kappa_var = rep(100000, d),
              prior_sigma_nu = d + nu_app, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
              prior_sigma_psi = diag(1, d)*(d+nu_app -d -1),

              # this section set the initial values of the parameters
              mu_init = mu_init,
              kappa_init = kappa_init,
              sigma_init = sigma_init,
              r_init =  r_init,
              

              # parameters for the adaptive part of Metropolis
              adapt_batch = batch_mcmc,
              adapt_a = a_mcmc,
              adapt_b = b_mcmc,
              adapt_alpha_target = alpha_target,
              sd_mu_scal = 1,
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
            kappa_out <- out_mcmc$kappa_out
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
              kappa_out[isim, ] <- kappa_out[isim, ] * diag(B)
              for (iobs in 1:n)
              {
                r_out[isim, iobs, ] <- r_out[isim, iobs, ] * diag(B)
              }
            }
# Bundle the true values, data, posterior samples, and diagnostics
# saved by the simulation study.
            res_list <- list(
              "runtime" = runtime,
              "n" = n,
              "d" = d,
              "kappa" = kappa,
              "mu" = mu,
              "Sigma_s" = Sigma_s,
              "Sigma_c" = Sigma_c,
              "x_c" = x_c,
              "x_s" = x_s,
              "theta" = theta,
              "mu_out" = mu_out,
              "kappa_out" = kappa_out,
              "sigma_s_out" = sigma_s_out,
              "sigma_c_out" = sigma_c_out,
              "pos_def_sigma" = out_mcmc$ess_acc_sigma
            )


            # Save one result file per scenario and chain.
            save(res_list, file = paste(
              "simulations/output/", name_sim, "tpn_simulations_results -",
              " select_d=", select_d,
              " select_n=", select_n,
              " select_kappa=", select_kappa,
              " select_sigma=", select_sigma,
              " select_chain=", select_chain,
              " seed_sigma=", select_sigma_type,
              " seed_data=" , select_data,
              ".Rdata"
            ))
            
            counter <- counter + 1
          }
        }
        
      }
    }
  }
}

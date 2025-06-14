rm(list = ls())

library(MCMCpack)
library(MASS)
library(truncnorm)
library(matrixcalc)
library(LaplacesDemon)
library(Rfast)
source("functions/general_functions.R")
source("functions/mcmc_tpn.R")


#### #### #### #### #### ####
#### Simulation
#### #### #### #### #### ####



# ========
# funciton to simulate a sigma which is valid after taking the absolute values
# ========
seed <- 2
out <- list()
parout <- list()
counter <- 1


for (select_n in 1:3)
{
  for (select_d in 1:3)
  {
    for (select_kappa in 1:4)
    {
      for (select_sigma in 1:2)
      {
        # Store the result in the list

        n <- c(20, 50, 100)[select_n] # number of observation
        d <- c(3, 6, 12)[select_d] # dimension of the torus
        dmax <- max(d)
        ### parameters
        mu <- c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6)[1:d]

        kappa_mat <- matrix(NA, ncol = dmax, nrow = 4)
        kappa_mat[1, ] <- rep(0.49, dmax)
        kappa_mat[2, ] <- rep(1.1, dmax)
        kappa_mat[3, ] <- rep(2.45, dmax)
        kappa_mat[4, ] <- rep(c(0.49, 1.1, 2.45), dmax / 3)
        kappa <- kappa_mat[select_kappa, 1:d]

        sigma_array <- array(0, c(dmax, dmax, 2))
        sigma_array[, , 1] <- diag(1, dmax)

        # 0.49   -> 0.3
        # 1,1     -> 0.6
        # 2.45   -> 0.9
        # simulation of sigma_s and sigma_c
        set.seed(1)
        d_test <- dmax
        max_attempts <- 100000
        attempts <- 0
        repeat{
          attempts <- attempts + 1
          Sigma_try <- sim_sigma(d_test + 1, diag(1, d_test))
          if (Sigma_try[[2]] || attempts >= max_attempts) {
            break
          }
        }
        sigma_array[, , 2] <- Sigma_try[[1]]
        Sigma_s <- sigma_array[1:d, 1:d, select_sigma]
        Sigma_c <- abs(Sigma_s)


        B <- diag(1 / diag(Sigma_s^0.5))
        Sigma_s <- B %*% Sigma_s %*% B
        Sigma_c <- B %*% Sigma_c %*% B


        set.seed(seed)

        # # # # # # # # # # # # # # # # # #
        # I simulate the linear variables
        # # # # # # # # # # # # # # # # # #
        x_c <- mvrnorm(n, kappa, Sigma_c)
        x_s <- mvrnorm(n, rep(0, d), Sigma_s)

        # # # # # # # # # # # # # # # # # #
        # And the circular ones
        # # # # # # # # # # # # # # # # # #
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

        ## Plot of the data
        ## This plot the marginal densities of the theta variables on the circle
        ## and the pairs of variables for d>1
        pdf(paste(
          "simulations/output/tpn_data - select_seed=", seed,
          " select_d=", select_d,
          " select_n=", select_n,
          " select_kappa=", select_kappa,
          " select_sigma=", select_sigma,
          ".pdf",
          sep = ""
        ))

        par(mfrow = c(2, 2))
        for (id in 1:d)
        {
          plot(density(theta[, id]),
            main = paste0("density of theta ", id),
            xlab = "theta", ylab = "density", xlim = c(0, 2 * pi)
          )
        }
        if (d > 1) {
          for (id in 1:(d - 1))
          {
            for (ij in (id + 1):d)
            {
              plot(theta[, id], theta[, ij],
                pch = 20, main = paste0("theta ", id, "vs theta ", ij),
                xlab = "theta", ylab = "theta", xlim = c(0, 2 * pi), ylim = c(0, 2 * pi)
              )
            }
          }
        }

        dev.off()

        #### #### #### #### #### ####
        #### MCMC function
        #### #### #### #### #### ####

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
          Sigma_c = Sigma_c,
          seed = seed
        )
        parout[[counter]] <- par_list
        save(parout, file = paste(
          "simulations/output/tpn_simulations_parameters -",
          " select_seed=", seed,
          # " select_d=", select_d,
          # " select_n=", select_n,
          # " select_kappa=", select_kappa,
          # " select_sigma=", select_sigma,
          ".Rdata"
        ))



        mmm <- 10
        out_mcmc <- mcmc_tpn(
          theta = theta, # the circualr data
          burnin = 1000 * mmm, # burnin
          thin = 1 * mmm, # thin
          iterations = 3000 * mmm, # total interations
          prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
          prior_mu_var = rep(100000, d),
          prior_kappa_mean = matrix(0, nrow = d, ncol = 1), # the prior for k is  TN(prior_kappa_mean,prior_kappa_var )
          prior_kappa_var = rep(100000, d),
          prior_sigma_nu = d + 2, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
          prior_sigma_psi = diag(1, d),

          # this section set the initial values of the parameters
          mu_init = rnorm(d) %% (2 * pi),
          kappa_init = abs(rnorm(d)),
          sigma_init = diag(1, d),
          r_init = matrix(1, nrow = n, ncol = d),

          # parameters for the adaptive part of Metropolis
          adapt_batch = 50,
          adapt_a = 1000,
          adapt_b = 1200,
          adapt_alpha_target = 0.234,
          sd_mu_scal = 1,
          par_sigma_adapt = 2000
        )

        # # # # # # # # # # # # # #
        # I extract the posterior samples of the parameters
        # # # # # # # # # # # # # #
        mu_out <- out_mcmc$mu_out
        kappa_out <- out_mcmc$kappa_out
        sigma_s_out <- out_mcmc$sigma_s_out
        sigma_c_out <- out_mcmc$sigma_c_out
        r_out <- out_mcmc$r_out
        ### identification
        mu_out <- mu_out %% (2 * pi)
        nsim <- nrow(mu_out)



        # # # # # # # # # # # # # #
        # The parameters must be indentified
        # # # # # # # # # # # # # #
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
        res_list <- list(
          "seed" = seed,
          "n" = n,
          "d" = d,
          "kappa" = kappa,
          "Sigma_s" = Sigma_s,
          "Sigma_c" = Sigma_c,
          "x_c" = x_c,
          "x_s" = x_s,
          "theta" = theta,
          "r" = r,
          "mu_out" = mu_out,
          "kappa_out" = kappa_out,
          "sigma_s_out" = sigma_s_out,
          "sigma_c_out" = sigma_c_out,
          "r_out" = r_out,
          "mcmc_sigma_s_out" = out_mcmc$sigma_s_out,
          "mcmc_sigma_c_out" = out_mcmc$sigma_c_out,
          "mcmc_r_out" = out_mcmc$r_out,
          "mcmc_kappa_out" = out_mcmc$kappa_out
        )



        out[[counter]] <- res_list


        save(out, file = paste(
          "simulations/output/tpn_simulations_results -",
          " select_seed=", seed,
          # " select_d=", select_d,
          # " select_n=", select_n,
          # " select_kappa=", select_kappa,
          # " select_sigma=", select_sigma,
          ".Rdata"
        ))


        ### plot of the parameters chain after identification with the true values
        pdf(paste(
          "simulations/output/tpn_chains - select_seed=", seed,
          " select_d=", select_d,
          " select_n=", select_n,
          " select_kappa=", select_kappa,
          " select_sigma=", select_sigma,
          ".pdf",
          sep = ""
        ))

        par(mfrow = c(3, 3))
        for (id in 1:d)
        {
          plot(mu_out[, id], type = "l", main = round(mu[id], 3))
          abline(h = mu[id], col = 2)
        }
        par(mfrow = c(3, 3))
        for (id in 1:d)
        {
          plot(kappa_out[, id], type = "l", main = round(kappa[id], 3))
          abline(h = kappa[id], col = 2)
        }
        par(mfrow = c(3, 3))
        h <- 1
        for (id in 1:d)
        {
          for (jd in 1:d)
          {
            plot(sigma_s_out[, h], type = "l", main = round(Sigma_s[id, jd], 3))
            abline(h = Sigma_s[id, jd], col = 2)
            h <- h + 1
          }
        }
        h <- 1
        for (id in 1:d)
        {
          for (jd in 1:d)
          {
            plot(sigma_c_out[, h], type = "l", main = round(Sigma_c[id, jd], 3))
            abline(h = Sigma_c[id, jd], col = 2)
            h <- h + 1
          }
        }
        dev.off()
        counter <- counter + 1
      }
    }
  }
}

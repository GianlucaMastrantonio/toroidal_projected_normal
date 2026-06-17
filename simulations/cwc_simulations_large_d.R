rm(list = ls())
library(CholWishart)

library(MCMCpack)
library(MASS)
library(truncnorm)
library(matrixcalc)
library(LaplacesDemon)
library(Rfast)
source("functions/general_functions.R")
source("functions/mcmc_cwn.R")
source("/beegfs/users/gmastrantonio/tokyo/codes/parameters_mcmc.R")
#### #### #### #### #### ####
#### Simulation
#### #### #### #### #### ####

# ========
# funciton to simulate a sigma which is valid after taking the absolute values
# ========


out <- list()
parout <- list()
counter <- 1

seed_list_chain <- c(1, 2)
seed_list_data <- 1:999
seed_list_sigma <- 1:999

args <- commandArgs(trailingOnly = TRUE)
app_d  = as.integer(args[1]) # 1:2
app_rho = as.integer(args[2]) # 1:4
app_sigma_ind_dep = as.integer(args[3]) # 1:2
app_chain = as.integer(args[4]) # 1:2
app_sigma = as.integer(args[5]) # 1:...
app_data = as.integer(args[6]) # 1:....
app_n = as.integer(args[7]) # 1:3
do_best_init = c(T,F)[as.integer(args[8])]

 #$d $k $sid $chain $sigma $data
name_sim <- paste(name_sim, "do_best_init=", do_best_init, sep = "")
select_d <- app_d

for (select_n in app_n:app_n)
{
  for (select_rho in app_rho:app_rho)
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
            # Store the result in the list

            d <- c(5, 50, 100)[select_d] # dimension of the torus
            n <- c(d*5, d*10, d*20)[select_n] # number of observation
            

            dmax <- max(d)
            ### parameters
            mu <- rep(0, d)
            kappa <- rep(0, d)


            sigma_array <- array(0, c(dmax, dmax, 2))
            sigma_array[, , 1] <- diag(1, dmax)


            # simulation of sigma_s and sigma_c
            set.seed(seed_sigma)
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

            set.seed(seed_data)

            # # # # # # # # # # # # # # # # # #
            # I simulate the linear variables
            # # # # # # # # # # # # # # # # # #
            x_c <- mvrnorm(n, kappa, Sigma_c)
            x_s <- mvrnorm(n, rep(0, d), Sigma_s)

            # # # # # # # # # # # # # # # # # #
            # And the circular ones
            # # # # # # # # # # # # # # # # # #
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

            ## Plot of the data
            ## This plot the marginal densities of the theta variables on the circle
            ## and the pairs of variables for d>1
            pdf(paste(
              "simulations/output/", name_sim, "cwc_data - ", 
              " select_d=", select_d,
              " select_n=", select_n,
              " select_rho=", select_rho,
              " select_sigma=", select_sigma,
              " select_chain=", select_chain,
              " seed_sigma=" ,select_sigma_type,
              " seed_data=", select_data,
              ".pdf",
              sep = ""
            ))

            par(mfrow = c(2, 2))
            for (id in 1:d)
            {
              plot(density(theta_cop[, id]),
                main = paste0("density of theta ", id),
                xlab = "theta", ylab = "density", xlim = c(0, 2 * pi)
              )
            }
            if (d > 1) {
              for (id in 1:(d - 1))
              {
                for (ij in (id + 1):d)
                {
                  plot(theta_cop[, id], theta_cop[, ij],
                    pch = 20, main = paste0("theta ", id, "vs theta ", ij),
                    xlab = "theta", ylab = "theta", xlim = c(0, 2 * pi), ylim = c(0, 2 * pi)
                  )
                }
              }
            }

            dev.off()
            ##### ##### ##### ##### ##### ##### #####
            ##### wrapped cauchy marginals
            ##### ##### ##### ##### ##### ##### #####


            mu <- rep(c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6), times = 20)[1:d]

            rho_mat <- matrix(NA, ncol = dmax, nrow = 4)
            rho_mat[1, ] <- rep(0.3, dmax) / 1
            rho_mat[2, ] <- rep(0.6, dmax) / 1
            rho_mat[3, ] <- rep(0.9, dmax) / 1
            rho_mat[4, ] <- (rep(c(0.3, 0.6, 0.9), (dmax + 3) / 3) / 1)[1:dmax]
            rho <- rho_mat[select_rho, 1:d]


            # cdf_wc = function(theta, mu, rho)
            # {
            #    n = length(theta)
            #    ret = rep(NA, n)
            #    d0 = cdf_wc_un(0, mu, rho)
            #    for(i in 1:n)
            #    {
            #        ret[i] = cdf_wc_un(theta[i], mu, rho)
            #    }
            #    return((ret-d0)%%1)
            # }


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
                # w <- which(cumulative_wrappedcauchy[[id]] > (theta_cop[iobs, id] / (2 * pi)))[1]
                ## theta[iobs, id] <- q_wc(theta_cop[iobs, id] / (2 * pi), 0, rho[id]) + mu[id]
                # theta[iobs, id] <- (theta_seq[w] + mu[id]) %% (2 * pi)
                theta[iobs, id] <- (q_wc(theta_cop[iobs, id] / (2 * pi), 0, rho[id]) + mu[id]) %% (2 * pi)
              }
            }

            pdf(paste(
              "simulations/output/", name_sim, "cwc_data - ", 
              " select_d=", select_d,
              " select_n=", select_n,
              " select_rho=", select_rho,
              " select_sigma=", select_sigma,
              " select_chain=", select_chain,
                " seed_sigma=" ,select_sigma_type,
              " seed_data=", select_data,
              ".pdf",
              sep = ""
            ))

            par(mfrow = c(2, 2))
            for (id in 1:d)
            {
              hist((theta[, id]),
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
            #parout[[counter]] <- par_list
            save(par_list, file = paste(
              "simulations/output/", name_sim, "cwc_simulations_parameters -",
              #" select_seed=", seed,
              " select_d=", select_d,
               " select_n=", select_n,
               " select_rho=", select_rho,
               " select_sigma=", select_sigma,
              " select_chain=", select_chain,
              " seed_sigma=" ,select_sigma_type,
              " seed_data=", select_data,
              ".Rdata"
            ))

            # rm(x_c)
            # rm(x_s)
            
            set.seed(seed_chain)
            
            
            if(do_best_init == TRUE)
            {
              mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x)), sum(cos(x))))
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
                  rrr[ik] <- sum(func_logd_wc(theta[, i], mu_init[i], rho_seq[ik]))
                }
                rho_init[i] <- rho_seq[which.max(rrr)]
                
                theta_app <- 2 * pi * func_cdf_wc(theta[, i] - mu_init[i], 0, rho_init[i])
                u <- runif(n)
                r_direct <- sqrt(-2 * log(u))
                r_app <- r_direct
                x_app <- r_app * cos(theta_app)
                y_app <- r_app * sin(theta_app)
                

                #ttt <- atan2(y_app, x_app)
                #qq <- quantile((ttt) + pi, prob = c(0.15, 0.85)) - pi
                
                #x_app <- x_app + kappa_init[i]
                #r_init[, i] <- sqrt(x_app^2 + y_app^2)
                #x_init[, i] <- x_app
                y_init[, i] <- y_app
              }
              sigma_init <- cov(y_init)
              sigma_init <- cov(y_init)
              while (inherits(try(chol(abs(sigma_init)), silent = TRUE), "try-error")) {

                sigma_init <- sigma_init + diag(0.01,d)

              }
            }else{
              mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x)), sum(cos(x))))
              rho_init <- runif(d, 0.5,0.9)
              r_init <- matrix(runif(n*d, 0.8,1.2), n, d)
              sigma_init <- diag(1, d)
            }
            mmm <- m_mcmc
            start <- Sys.time()
            out_mcmc <- mcmc_cwc(
              theta = theta, # the circualr data
              burnin = burnin_mcmc * mmm, # burnin
              thin = thin_mcmc * mmm, # thin
              iterations = iter_mcmc * mmm, # total interations
              #burnin =10, # burnin
              #thin = 1 , # thin
              #iterations = 100, # total interations
              prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
              prior_mu_var = rep(100000, d),
              prior_rho_a = rep(1, d), # the prior for B()
              prior_rho_b = rep(1, d),
              prior_sigma_nu = d + 2, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
              prior_sigma_psi = diag(1, d),


              # this section set the initial values of the parameters
              mu_init = mu_init,
              rho_init = rho_init,
              sigma_init = sigma_init,
              r_init =r_init,

              # parameters for the adaptive part of Metropolis
                adapt_batch = batch_mcmc,
              adapt_a = a_mcmc,
              adapt_b = b_mcmc,
              adapt_alpha_target = alpha_target,
              sd_mu_scal = 1,
              sd_rho_scal = 0.1,
              par_sigma_adapt = par_adapt_mcmc
            )
            end <- Sys.time()
            runtime <- end - start
            # # # # # # # # # # # # # #
            # I extract the posterior samples of the parameters
            # # # # # # # # # # # # # #
            mu_out <- out_mcmc$mu_out
            rho_out <- out_mcmc$rho_out
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
              # kappa_out[isim, ] <- kappa_out[isim, ] * diag(B)
              # for (iobs in 1:n)
              # {
              #  r_out[isim, iobs, ] <- r_out[isim, iobs, ] * diag(B)
              # }
            }
            res_list <- list(
              #"wc seed" = seed,
              "runtime" = runtime,
              "n" = n,
              "d" = d,
              "rho" = rho,
              "Sigma_s" = Sigma_s,
              "Sigma_c" = Sigma_c,
              "x_c" = x_c,
              "x_s" = x_s,
              "theta" = theta,
              "r" = r,
              "mu_out" = mu_out,
              "sigma_s_out" = sigma_s_out,
              "sigma_c_out" = sigma_c_out,
              "r_out" = r_out,
              "mcmc_sigma_s_out" = out_mcmc$sigma_s_out,
              "mcmc_sigma_c_out" = out_mcmc$sigma_c_out,
              "rho_out" = rho_out,
              "pos_def_sigma" = out_mcmc$save_acc_sigma
            )
            #out[[counter]] <- res_list

            save(res_list, file = paste(
              "simulations/output/", name_sim, "cwc_simulations_results -",
              #" select_seed=", seed,
              " select_d=", select_d,
               " select_n=", select_n,
               " select_rho=", select_rho,
               " select_sigma=", select_sigma,
              " select_chain=", select_chain,
                " seed_sigma=" ,select_sigma_type,
              " seed_data=", select_data,
              ".Rdata"
            ))


            ### plot of the parameters chain after identification with the true values
            pdf(paste(
              "simulations/output/", name_sim, "cwc_chains -  ",  
              " select_d=", select_d,
              " select_n=", select_n,
              " select_rho=", select_rho,
              " select_sigma=", select_sigma,
              " select_chain=", select_chain,
              " seed_sigma=" ,select_sigma_type,
              " seed_data=", select_data,
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
              plot(rho_out[, id], type = "l", main = round(rho[id], 3))
              abline(h = rho[id], col = 2)
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
  }
}

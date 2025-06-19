mcmc_cwc <- function(
    theta,
    burnin,
    thin,
    iterations,
    prior_mu_mean,
    prior_mu_var,
    prior_rho_a,
    prior_rho_b,
    prior_sigma_nu,
    prior_sigma_psi,
    mu_init,
    rho_init,
    sigma_init,
    r_init,
    adapt_batch,
    adapt_a,
    adapt_b,
    adapt_alpha_target,
    sd_mu_scal,
    sd_rho_scal,
    par_sigma_adapt,
    na_index = list(NA)) {
    d <- dim(theta)[2]
    n <- dim(theta)[1]
    sample_to_save <- round((iterations - burnin) / thin)

    # This containt the posterior sample that we are going to save
    mu_out <- matrix(NA, nrow = sample_to_save, ncol = d)
    rho_out <- matrix(NA, nrow = sample_to_save, ncol = d)
    sigma_s_out <- matrix(NA, nrow = sample_to_save, ncol = d^2)
    sigma_c_out <- matrix(NA, nrow = sample_to_save, ncol = d^2)
    r_out <- array(NA, c(sample_to_save, n, d))

    # these objects contains the corrent value of the parameters
    mu_mcmc <- matrix(NA, ncol = 1, nrow = d)
    rho_mcmc <- matrix(NA, ncol = 1, nrow = d)
    sigma_s_mcmc <- matrix(NA, nrow = d, ncol = d)
    sigma_c_mcmc <- matrix(NA, nrow = d, ncol = d)
    lambda_c_mcmc <- matrix(NA, nrow = d, ncol = d)
    lambda_s_mcmc <- matrix(NA, nrow = d, ncol = d)
    r_mcmc <- matrix(abs(rnorm(n * d)), nrow = n, ncol = d)
    # for(iobs in 1:n)
    # {
    #    for(id in 1:d)
    #    {
    #        r_mcmc[iobs, id] = (x_c[iobs, id]^2+x_s[iobs, id]^2)^0.5
    #    }
    # }

    # the linear variables
    x_s_mcmc <- matrix(1, nrow = n, ncol = d)
    x_c_mcmc <- matrix(1, nrow = n, ncol = d)

    mu_mcmc[] <- mu_init
    rho_mcmc[] <- rho_init
    r_mcmc[] <- r_init



    # the adaptive part of the Metropolis
    sd_r <- matrix(1, nrow = n, ncol = d)
    alpha_r <- matrix(0, nrow = n, ncol = d)
    sd_mu <- matrix(1, nrow = d, ncol = 1) * sd_mu_scal
    alpha_mu <- matrix(0, nrow = d, ncol = 1)
    sd_rho <- matrix(1, nrow = d, ncol = 1) * sd_rho_scal
    alpha_rho <- matrix(0, nrow = d, ncol = 1)
    alpha_sigma <- 0




    # other objects that containt the current value
    sigma_s_mcmc[, ] <- (sigma_init + t(sigma_init)) / 2 # i did this because soimethimes, the matrices are not exactly simmetrical
    sigma_c_mcmc[, ] <- abs(sigma_s_mcmc)

    lambda_c_mcmc <- solve(sigma_c_mcmc)
    lambda_s_mcmc <- solve(sigma_s_mcmc)

    cond_mean_s <- matrix(NA, nrow = n, ncol = d)
    cond_mean_c <- matrix(NA, nrow = n, ncol = d)
    # ====
    # Missings
    # ====
    there_are_na <- !is.na(na_index[[1]][[1]])
    missig_out <- list()
    if (there_are_na) {
        for (id in 1:d)
        {
            missig_out[[id]] <- matrix(NA, nrow = sample_to_save, ncol = length(na_index[[id]]))
        }
    } else {
        for (id in 1:d)
        {
            missig_out[[id]] <- matrix(NA, nrow = sample_to_save, ncol = 1)
        }
    }

    tf_missig_out <- matrix(T, nrow = n, ncol = d)

    theta_cop <- theta
    theta_cop_prop <- theta

    for (id in 1:d)
    {
        for (iobs in 1:n)
        {
            theta_cop[iobs, id] <- 2 * pi * func_cdf_wc(theta[iobs, id] - mu_mcmc[id], 0, rho_mcmc[id])
        }
    }
    for (id in 1:d)
    {
        x_s_mcmc[, id] <- r_mcmc[, id] * sin(theta_cop[, id] - 0)
        x_c_mcmc[, id] <- r_mcmc[, id] * cos(theta_cop[, id] - 0)
    }
    x_s_prop <- x_s_mcmc
    x_c_prop <- x_c_mcmc

    x_s_zeta <- matrix(NA, nrow = d)
    x_c_zeta <- matrix(NA, nrow = d)
    theta_zeta <- matrix(NA, nrow = d)
    sum_iter <- 0
    burn_thin <- burnin
    # * WAIC
    delta_copulat <- seq(0, 2 * pi, length.out = 10 * 361)[2]

    vec_zero <- matrix(0, nrow = d)
    sum_log_dens_data <- rep(0, n)
    sum_dens_data <- rep(0, n)

    n_samp_cop <- 10000
    samp_wc <- matrix(NA, nrow = n_samp_cop, ncol = d)
    for (imcmc in 1:sample_to_save)
    {
        for (jmcmc in 1:burn_thin)
        {
            ###### all the step follow what I wrote on the tex file


            sum_iter <- sum_iter + 1
            ## obs-specific parameters
            if ((sum_iter %% 50) == 0) {
                print(paste("Iteration:", sum_iter))
            }
            ### Inversi of PN varianbles
            # theta_seq <- seq(0, 2 * pi, by = delta_copulat)
            # cumulative_wrappedcauchy <- array(NA, c(d, length(theta_seq)))
            # for (id in 1:d)
            # {
            #    cumulative_wrappedcauchy[id, ] <- func_cdf_wc(theta_seq, mu_mcmc[id], rho_mcmc[id])
            # }
            ## missing
            if (there_are_na) {
                x_s_prop <- x_s_mcmc
                x_c_prop <- x_c_mcmc
                theta_cop_prop <- theta_cop
                for (id in 1:d)
                {
                    for (iobs in na_index[[id]])
                    {
                        cond_var_c <- 1 / lambda_c_mcmc[id, id]
                        cond_var_s <- 1 / lambda_s_mcmc[id, id]

                        cond_mean_c <- 0 - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - 0))
                        cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])


                        theta_prop <- rnorm(1, theta[iobs, id], sample(c(0.5, 0.1, 0.05, 0.01), 1)) %% (2 * pi)


                        theta_cop_prop[iobs, id] <- 2 * pi * func_cdf_wc(theta_prop - mu_mcmc[id], 0, rho_mcmc[id])
                        x_s_prop[iobs, id] <- r_mcmc[iobs, id] * sin(theta_cop_prop[iobs, id])
                        x_c_prop[iobs, id] <- r_mcmc[iobs, id] * cos(theta_cop_prop[iobs, id])

                        mh_ratio <- 0

                        mh_ratio <- mh_ratio + dnorm(x_c_prop[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                        mh_ratio <- mh_ratio + dnorm(x_s_prop[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

                        mh_ratio <- mh_ratio - dnorm(x_c_mcmc[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                        mh_ratio <- mh_ratio - dnorm(x_s_mcmc[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

                        mh_ratio <- mh_ratio + log(func_d_wc(theta_prop, mu_mcmc[id], rho_mcmc[id]))

                        mh_ratio <- mh_ratio - log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]))

                        if (is.na(mh_ratio)) {
                            print("mh_ratio is NA missing")
                            mh_ratio <- log(0)
                        }

                        if (runif(1, 0, 1) < exp(mh_ratio)) {
                            theta[iobs, id] <- theta_prop
                            x_s_mcmc[iobs, id] <- x_s_prop[iobs, id]
                            x_c_mcmc[iobs, id] <- x_c_prop[iobs, id]
                            theta_cop[iobs, id] <- theta_cop_prop[iobs, id]
                        } else {
                            x_s_prop[iobs, id] <- x_s_mcmc[iobs, id]
                            x_c_prop[iobs, id] <- x_c_mcmc[iobs, id]
                            theta_cop_prop[iobs, id] <- theta_cop[iobs, id]
                        }
                    }
                }
            }
            for (id in 1:d)
            {
                # print(theta_zeta[id, k])
                # print(theta[iobs, id])
                # print(mu_mcmc[id, k])
                # print(rho_mcmc[id, k])
                theta_zeta[id] <- 2 * pi * func_cdf_wc(theta[iobs, id] - mu_mcmc[id], 0, rho_mcmc[id])
                x_s_zeta[id] <- r_mcmc[iobs, id] * sin(theta_zeta[id] - 0)
                x_c_zeta[id] <- r_mcmc[iobs, id] * cos(theta_zeta[id] - 0)
            }
            #### mu
            x_s_prop <- x_s_mcmc
            x_c_prop <- x_c_mcmc
            theta_cop_prop <- theta_cop
            for (id in 1:d)
            {
                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]

                mu_prop <- rnorm(1, mu_mcmc[id], sd_mu[id])



                for (iobs in 1:n)
                {
                    theta_cop_prop[iobs, id] <- 2 * pi * func_cdf_wc(theta[iobs, id] - mu_prop, 0, rho_mcmc[id])
                }
                x_s_prop[, id] <- r_mcmc[, id] * sin(theta_cop_prop[, id])
                x_c_prop[, id] <- r_mcmc[, id] * cos(theta_cop_prop[, id])

                mh_ratio <- 0
                for (iobs in 1:n)
                {
                    cond_mean_c <- 0 - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - 0))
                    cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

                    mh_ratio <- mh_ratio + dnorm(x_c_prop[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                    mh_ratio <- mh_ratio + dnorm(x_s_prop[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

                    mh_ratio <- mh_ratio - dnorm(x_c_mcmc[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                    mh_ratio <- mh_ratio - dnorm(x_s_mcmc[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

                    mh_ratio <- mh_ratio + log(func_d_wc(theta[iobs, id], mu_prop, rho_mcmc[id]))

                    mh_ratio <- mh_ratio - log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]))
                }

                mh_ratio <- mh_ratio + dnorm(mu_prop, prior_mu_mean[id], prior_mu_var[id]^0.5, log = T)
                mh_ratio <- mh_ratio - dnorm(mu_mcmc[id], prior_mu_mean[id], prior_mu_var[id]^0.5, log = T)

                alpha_mh <- min(1, exp(mh_ratio))
                alpha_mu[id] <- alpha_mu[id] + alpha_mh

                # if(id == 1)
                # {
                #    print(alpha_mh)
                #    print(sd_mu[id])
                # }
                if (is.na(alpha_mh)) {
                    print("mh_ratio is NA mu")
                    alpha_mh <- log(0)
                }
                if (runif(1, 0, 1) < alpha_mh) {
                    mu_mcmc[id] <- mu_prop
                    x_s_mcmc[, id] <- x_s_prop[, id]
                    x_c_mcmc[, id] <- x_c_prop[, id]
                    theta_cop[, id] <- theta_cop_prop[, id]
                } else {
                    x_s_prop[, id] <- x_s_mcmc[, id]
                    x_c_prop[, id] <- x_c_mcmc[, id]
                    theta_cop_prop[, id] <- theta_cop[, id]
                }
            }
            #### rho
            x_s_prop <- x_s_mcmc
            x_c_prop <- x_c_mcmc
            theta_cop_prop <- theta_cop
            for (id in 1:d)
            {
                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]

                rho_mcmc_app <- log(rho_mcmc[id] / (1 - rho_mcmc[id]))
                rho_prop_app <- rnorm(1, rho_mcmc_app, sd_rho[id])
                # rho_prop_app = rho_mcmc_app
                rho_prop <- exp(rho_prop_app) / (1 + exp(rho_prop_app))


                for (iobs in 1:n)
                {
                    theta_cop_prop[iobs, id] <- 2 * pi * func_cdf_wc(theta[iobs, id] - mu_mcmc[id], 0, rho_prop)
                }
                x_s_prop[, id] <- r_mcmc[, id] * sin(theta_cop_prop[, id])
                x_c_prop[, id] <- r_mcmc[, id] * cos(theta_cop_prop[, id])

                mh_ratio <- 0
                for (iobs in 1:n)
                {
                    cond_mean_c <- 0 - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - 0))
                    cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

                    mh_ratio <- mh_ratio + dnorm(x_c_prop[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                    mh_ratio <- mh_ratio + dnorm(x_s_prop[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

                    mh_ratio <- mh_ratio - dnorm(x_c_mcmc[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                    mh_ratio <- mh_ratio - dnorm(x_s_mcmc[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

                    mh_ratio <- mh_ratio + log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_prop))

                    mh_ratio <- mh_ratio - log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]))
                }

                mh_ratio <- mh_ratio + (dbeta(rho_prop, prior_rho_a[id], prior_rho_b[id], log = T) + rho_prop_app - 2 * log(1 + exp(rho_prop_app)))
                mh_ratio <- mh_ratio - (dbeta(rho_mcmc[id], prior_rho_a[id], prior_rho_b[id], log = T) + rho_mcmc_app - 2 * log(1 + exp(rho_mcmc_app)))
                # print(round(mh_ratio,4))

                alpha_mh <- min(1, exp(mh_ratio))
                alpha_rho[id] <- alpha_rho[id] + alpha_mh


                # rho = exp(x)/(1+exp(x))

                # f_x(x) = f_r(r) dr/dx
                if (is.na(alpha_mh)) {
                    print("mh_ratio is NA rho")
                    alpha_mh <- log(0)
                }
                if (runif(1, 0, 1) < alpha_mh) {
                    rho_mcmc[id] <- rho_prop
                    x_s_mcmc[, id] <- x_s_prop[, id]
                    x_c_mcmc[, id] <- x_c_prop[, id]
                    theta_cop[, id] <- theta_cop_prop[, id]
                } else {
                    x_s_prop[, id] <- x_s_mcmc[, id]
                    x_c_prop[, id] <- x_c_mcmc[, id]
                    theta_cop_prop[, id] <- theta_cop[, id]
                }
            }








            ##### r
            for (id in 1:d)
            {
                x_s_prop <- x_s_mcmc
                x_c_prop <- x_c_mcmc

                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]
                cond_sigma <- diag(c(1 / lambda_c_mcmc[id, id], 1 / lambda_s_mcmc[id, id]))
                inv_cond_sigma <- solve(cond_sigma)
                for (iobs in 1:n)
                {
                    cond_mean_c <- 0 - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - 0))
                    cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

                    cond_mean <- matrix(c(cond_mean_c, cond_mean_s), nrow = 2)

                    uvec <- matrix(c(cos(theta_cop[iobs, id] - 0), sin(theta_cop[iobs, id] - 0)), nrow = 2)


                    A <- t(uvec) %*% (inv_cond_sigma) %*% uvec
                    B <- t(uvec) %*% (inv_cond_sigma) %*% cond_mean


                    v1sim <- runif(1, 0, exp(-0.5 * A * (r_mcmc[iobs, id] - B / A)^2))
                    v2sim <- runif(1, 0, 1)
                    rho1 <- B / A + max(-B / A, -sqrt(-2 * log(v1sim) / A))
                    rho2 <- B / A + sqrt(-2 * log(v1sim) / A)

                    r_mcmc[iobs, id] <- ((rho2^2 - rho1^2) * v2sim + rho1^2)^(1 / 2)

                    x_s_mcmc[iobs, id] <- r_mcmc[iobs, id] * sin(theta_cop[iobs, id] - 0)
                    x_c_mcmc[iobs, id] <- r_mcmc[iobs, id] * cos(theta_cop[iobs, id] - 0)
                }
            }

            ### sigma
            # print(cov(x_s_mcmc[, 1], x_s_mcmc[, 3]))
            # print(c(x_s_mcmc[1, 1:6]))
            # print(theta_cop[1, ])
            # print(r_mcmc[1, 1:6])


            # error("")
            nu <- par_sigma_adapt + d + 1
            prop_sigma <- rinvwishart(nu, sigma_s_mcmc * (nu - (d + 1)))

            test_sigma(prop_sigma)

            # if ((test_sigma(prop_sigma) == TRUE) & (test_sigma2(prop_sigma) == TRUE)) {
            if (test_sigma(prop_sigma) == TRUE) {
                sigma_s_prop <- prop_sigma
                sigma_c_prop <- abs(prop_sigma)

                lambda_c_prop <- solve(sigma_c_prop)
                lambda_s_prop <- solve(sigma_s_prop)

                log_det_c_mcmc <- determinant(sigma_c_mcmc, logarithm = T)
                log_det_c_prop <- determinant(sigma_c_prop, logarithm = T)

                log_det_s_mcmc <- determinant(sigma_s_mcmc, logarithm = T)
                log_det_s_prop <- determinant(sigma_s_prop, logarithm = T)


                mh_ratio <- 0
                for (iobs in 1:n)
                {
                    mh_ratio <- mh_ratio + (-0.5 * c(log_det_c_prop$modulus) - 0.5 * t(x_c_mcmc[iobs, ] - 0) %*% lambda_c_prop %*% (x_c_mcmc[iobs, ] - 0))
                    mh_ratio <- mh_ratio + (-0.5 * c(log_det_s_prop$modulus) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_prop %*% (x_s_mcmc[iobs, ]))

                    mh_ratio <- mh_ratio - (-0.5 * c(log_det_c_mcmc$modulus) - 0.5 * t(x_c_mcmc[iobs, ] - 0) %*% lambda_c_mcmc %*% (x_c_mcmc[iobs, ] - 0))
                    mh_ratio <- mh_ratio - (-0.5 * c(log_det_s_mcmc$modulus) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_mcmc %*% (x_s_mcmc[iobs, ]))
                }
                # prior
                mh_ratio <- mh_ratio + dinvwishart(sigma_s_prop, prior_sigma_nu, prior_sigma_psi, log = T)
                mh_ratio <- mh_ratio - dinvwishart(sigma_s_mcmc, prior_sigma_nu, prior_sigma_psi, log = T)

                ## proposal
                mh_ratio <- mh_ratio - dinvwishart(sigma_s_prop, nu, sigma_s_mcmc * (nu - (d + 1)), log = T)
                mh_ratio <- mh_ratio + dinvwishart(sigma_s_mcmc, nu, sigma_s_prop * (nu - (d + 1)), log = T)

                alpha_sigma <- alpha_sigma + min(1, exp(mh_ratio))

                if (is.na(mh_ratio)) {
                    print("mh_ratio is NA")
                    mh_ratio <- log(0)
                }
                if (runif(1, 0, 1) < exp(mh_ratio)) {
                    # print("ACC")
                    sigma_s_mcmc <- sigma_s_prop
                    sigma_c_mcmc <- sigma_c_prop

                    lambda_s_mcmc <- solve(sigma_s_mcmc)
                    lambda_c_mcmc <- solve(sigma_c_mcmc)
                }
            }




            ### update of the adaptive parameters
            if ((sum_iter %% adapt_batch == 0)) {
                alpha_mu <- alpha_mu / adapt_batch
                alpha_rho <- alpha_rho / adapt_batch
                # alpha_r = alpha_r/adapt_batch
                alpha_sigma <- alpha_sigma / adapt_batch
                # print(cbind(alpha_rho,sd_rho))
                # print(par_sigma_adapt)
                if ((sum_iter < burnin)) {
                    for (id in 1:d)
                    {
                        sd_mu[id] <- exp(log(sd_mu[id]) + adapt_a / (adapt_b + sum_iter) * (alpha_mu[id] - adapt_alpha_target))
                        alpha_mu[id] <- 0

                        sd_rho[id] <- exp(log(sd_rho[id]) + adapt_a / (adapt_b + sum_iter) * (alpha_rho[id] - adapt_alpha_target))
                        alpha_rho[id] <- 0
                        # for(iobs in 1:n)
                        # {
                        #    sd_r[iobs,id] = exp(log(sd_r[iobs,id]) +  adapt_a/(adapt_b+sum_iter)*(alpha_r[iobs,id] - adapt_alpha_target) )
                        #    alpha_r[iobs,id] = 0
                        # }
                    }

                    # sigma
                    par_sigma_adapt <- exp(log(par_sigma_adapt) - adapt_a / (adapt_b + sum_iter) * (alpha_sigma - adapt_alpha_target))
                    alpha_sigma <- 0
                }
            }
        }
        burn_thin <- thin

        # i save the current values of the chains
        mu_out[imcmc, ] <- mu_mcmc
        rho_out[imcmc, ] <- rho_mcmc
        sigma_s_out[imcmc, ] <- c(sigma_s_mcmc)
        sigma_c_out[imcmc, ] <- c(sigma_c_mcmc)
        r_out[imcmc, , ] <- r_mcmc

        if (there_are_na) {
            for (id in 1:d)
            {
                if (length(na_index[[id]]) > 0) {
                    missig_out[[id]][imcmc, ] <- theta[na_index[[id]], id]
                }
            }
        }
        for (iobs in 1:n)
        {
            app <- dmvnorm(x_c_zeta[tf_missig_out[iobs, ]], vec_zero[tf_missig_out[iobs, ]], sigma_c_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)

            app <- app + dmvnorm(x_s_zeta[tf_missig_out[iobs, ]], vec_zero[tf_missig_out[iobs, ]], sigma_s_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)

            for (id in 1:d)
            {
                if (tf_missig_out[iobs, id] == T) {
                    app <- app + func_logd_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]) + log(r_mcmc[iobs, id]) + log(2 * pi)
                }
            }


            sum_log_dens_data[iobs] <- sum_log_dens_data[iobs] + app
            sum_dens_data[iobs] <- sum_dens_data[iobs] + exp(app)
        }
    }
    waic_llpd <- sum(log(sum_dens_data[iobs] / sample_to_save))

    p_waic <- 2 * sum(log(sum_dens_data[iobs] / sample_to_save) - sum_log_dens_data[iobs] / sample_to_save)
    return(list(mu_out = mu_out, rho_out = rho_out, sigma_s_out = sigma_s_out, sigma_c_out = sigma_c_out, r_out = r_out, missig_out = missig_out, waic = 2 * (waic_llpd - p_waic)))
}

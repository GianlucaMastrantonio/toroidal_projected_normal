mcmc_tpn <- function(
    theta,
    burnin,
    thin,
    iterations,
    prior_mu_mean,
    prior_mu_var,
    prior_kappa_mean,
    prior_kappa_var,
    prior_sigma_nu,
    prior_sigma_psi,
    mu_init,
    kappa_init,
    sigma_init,
    r_init,
    adapt_batch,
    adapt_a,
    adapt_b,
    adapt_alpha_target,
    sd_mu_scal,
    par_sigma_adapt,
    na_index = list(NA)) {
    d <- dim(theta)[2]
    n <- dim(theta)[1]
    sample_to_save <- round((iterations - burnin) / thin)

    # This containt the posterior sample that we are going to save
    mu_out <- matrix(NA, nrow = sample_to_save, ncol = d)
    kappa_out <- matrix(NA, nrow = sample_to_save, ncol = d)
    sigma_s_out <- matrix(NA, nrow = sample_to_save, ncol = d^2)
    sigma_c_out <- matrix(NA, nrow = sample_to_save, ncol = d^2)
    r_out <- array(NA, c(sample_to_save, n, d))

    # these objects contains the corrent value of the parameters
    mu_mcmc <- matrix(NA, ncol = 1, nrow = d)
    kappa_mcmc <- matrix(NA, ncol = 1, nrow = d)
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
    kappa_mcmc[] <- kappa_init
    r_mcmc[] <- r_init

    for (id in 1:d)
    {
        x_s_mcmc[, id] <- r_mcmc[, id] * sin(theta[, id] - mu_mcmc[id])
        x_c_mcmc[, id] <- r_mcmc[, id] * cos(theta[, id] - mu_mcmc[id])
    }
    x_s_prop <- x_s_mcmc
    x_c_prop <- x_c_mcmc

    # the adaptive part of the Metropolis
    sd_r <- matrix(1, nrow = n, ncol = d)
    alpha_r <- matrix(0, nrow = n, ncol = d)
    sd_mu <- matrix(1, nrow = d, ncol = 1) * sd_mu_scal
    alpha_mu <- matrix(0, nrow = d, ncol = 1)
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
    # par_1 = n+prior_sigma_nu
    # par_2 = prior_sigma_psi
    # for(iobs in 1:n)
    # {
    #    par_2 = par_2 + t(x_s_mcmc)%*%(x_s_mcmc)
    # }
    # max_attempts = 1000
    # attempts = 0
    # repeat
    # {
    #    attempts = attempts + 1
    #    Sigma_try = sim_sigma(par_1,par_2)
    #    if (Sigma_try[[2]]) {
    #      break
    #    }
    #    if(attempts >= max_attempts)
    #    {
    #        error("too many iterations")
    #    }
    # }

    # sigma_s_mcmc[,] = Sigma_try[[1]]
    # sigma_c_mcmc[,] = abs(sigma_s_mcmc)

    # lambda_c_mcmc = solve(sigma_c_mcmc)
    # lambda_s_mcmc = solve(sigma_c_mcmc)
    sum_iter <- 0
    burn_thin <- burnin
    # * WAIC
    vec_zero <- matrix(0, nrow = d)
    sum_log_dens_data <- rep(0, n)
    sum_dens_data <- rep(0, n)
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

            #### missing
            if (there_are_na) {
                for (id in 1:d)
                {
                    for (iobs in na_index[[id]])
                    {
                        cond_var_c <- 1 / lambda_c_mcmc[id, id]
                        cond_var_s <- 1 / lambda_s_mcmc[id, id]

                        cond_mean_c <- kappa_mcmc[id] - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - kappa_mcmc[-id]))
                        cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

                        x_c_mcmc[iobs, id] <- rnorm(1, cond_mean_c, cond_var_c^0.5)
                        x_s_mcmc[iobs, id] <- rnorm(1, cond_mean_s, cond_var_s^0.5)

                        r_mcmc[iobs, id] <- sqrt((x_c_mcmc[iobs, id]^2 + x_s_mcmc[iobs, id]^2))
                        theta[iobs, id] <- atan2(x_s_mcmc[iobs, id], x_c_mcmc[iobs, id]) + mu_mcmc[id]
                    }
                }
            }

            ### mu
            x_s_prop <- x_s_mcmc
            x_c_prop <- x_c_mcmc
            for (id in 1:d)
            {
                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]

                mu_prop <- rnorm(1, mu_mcmc[id], sd_mu[id])

                x_s_prop[, id] <- r_mcmc[, id] * sin(theta[, id] - mu_prop)
                x_c_prop[, id] <- r_mcmc[, id] * cos(theta[, id] - mu_prop)

                mh_ratio <- 0
                for (iobs in 1:n)
                {
                    cond_mean_c <- kappa_mcmc[id] - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - kappa_mcmc[-id]))
                    cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

                    mh_ratio <- mh_ratio + dnorm(x_c_prop[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                    mh_ratio <- mh_ratio + dnorm(x_s_prop[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

                    mh_ratio <- mh_ratio - dnorm(x_c_mcmc[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
                    mh_ratio <- mh_ratio - dnorm(x_s_mcmc[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)
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
                if (runif(1, 0, 1) < alpha_mh) {
                    mu_mcmc[id] <- mu_prop
                    x_s_mcmc[, id] <- x_s_prop[, id]
                    x_c_mcmc[, id] <- x_c_prop[, id]
                } else {
                    x_s_prop[, id] <- x_s_mcmc[, id]
                    x_c_prop[, id] <- x_c_mcmc[, id]
                }
            }




            #### kappa
            for (id in 1:d)
            {
                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]

                sum_app_x <- 0
                for (iobs in 1:n)
                {
                    sum_app_x <- sum_app_x + x_c_mcmc[iobs, id] - (-cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - kappa_mcmc[-id])))
                }


                par_2 <- 1 / (n / cond_var_c + 1 / prior_kappa_var[id])
                par_1 <- par_2 * (sum_app_x / cond_var_c + prior_kappa_mean[id] / prior_kappa_var[id])

                kappa_mcmc[id] <- rtruncnorm(1, a = 0, b = Inf, mean = par_1, sd = par_2^0.5)
            }


            ## mu and k
            mu_prop <- rep(NA, d)
            kappa_prop <- rep(NA, d)

            mu_p_c <- matrix(0, nrow = d)
            var_p_c <- matrix(0, nrow = d, ncol = d)
            mu_p_s <- matrix(0, nrow = d)
            var_p_s <- matrix(0, nrow = d, ncol = d)
            # for (iobs in 1:n)
            # {
            #    x_c_prop[iobs, ] <- r_mcmc[iobs, ] * cos(theta[iobs, ])
            #    x_s_prop[iobs, ] <- r_mcmc[iobs, ] * sin(theta[iobs, ])
            # }
            x_c_prop[, ] <- r_mcmc * cos(theta)
            x_s_prop[, ] <- r_mcmc * sin(theta)


            var_p_c <- solve(n * lambda_c_mcmc)
            var_p_s <- solve(n * lambda_s_mcmc)
            for (iobs in 1:n)
            {
                mu_p_c <- mu_p_c + lambda_c_mcmc %*% t(x_c_prop[iobs, , drop = F])
                mu_p_s <- mu_p_s + lambda_s_mcmc %*% t(x_s_prop[iobs, , drop = F])
            }
            mu_p_c <- var_p_c %*% mu_p_c
            mu_p_s <- var_p_s %*% mu_p_s

            prop_s <- mu_p_s + t(chol(var_p_s)) %*% rnorm(d)
            prop_c <- mu_p_c + t(chol(var_p_c)) %*% rnorm(d)

            for (id in 1:d)
            {
                mu_prop[id] <- atan2(prop_s[id], prop_c[id])
                kappa_prop[id] <- (prop_s[id]^2 + prop_c[id]^2)^0.5
            }


            mh_ratio <- sum(log(kappa_prop)) - sum(log(kappa_mcmc))

            mh_ratio <- mh_ratio + sum(dnorm(kappa_prop, prior_kappa_mean, prior_kappa_var^0.5, log = T)) - sum(dnorm(kappa_mcmc, prior_kappa_mean, prior_kappa_var^0.5, log = T))

            if (log(runif(1, 0, 1)) < (mh_ratio)) {
                mu_mcmc <- mu_prop
                kappa_mcmc <- kappa_prop
                for (iobs in 1:n)
                {
                    x_c_mcmc[iobs, ] <- r_mcmc[iobs, ] * cos(theta[iobs, ] - mu_mcmc)
                    x_s_mcmc[iobs, ] <- r_mcmc[iobs, ] * sin(theta[iobs, ] - mu_mcmc)
                }
            }
            ##### r

            # for(id in 1:d)
            # {


            #    x_s_prop = x_s_mcmc
            #    x_c_prop = x_c_mcmc

            #    cond_var_c = 1/lambda_c_mcmc[id,id]
            #    cond_var_s = 1/lambda_s_mcmc[id,id]
            #    for(iobs in 1:n)
            #    {
            #        cond_mean_c = kappa_mcmc[id] - cond_var_c*sum(lambda_c_mcmc[id,-id]*(x_c_mcmc[iobs,-id] - kappa_mcmc[-id]))
            #        cond_mean_s =  - cond_var_s*sum(lambda_s_mcmc[id,-id]*x_s_mcmc[iobs,-id])

            #        r_prop = exp(rnorm(1,log(r_mcmc[iobs, id]) ,  sd_r[iobs, id] ))

            #        x_s_prop[iobs,id] = r_prop*sin(theta[iobs,id] - mu_mcmc[id])
            #        x_c_prop[iobs,id] = r_prop*cos(theta[iobs,id] - mu_mcmc[id])

            #        mh_ratio = 0


            #        mh_ratio = mh_ratio + dnorm(x_s_prop[iobs,id], cond_mean_s, cond_var_s^0.5, log=T)
            #        mh_ratio = mh_ratio + dnorm(x_c_prop[iobs,id], cond_mean_c, cond_var_c^0.5, log=T)
            #        mh_ratio = mh_ratio + 2*log(r_prop)

            #        mh_ratio = mh_ratio - dnorm(x_s_mcmc[iobs,id], cond_mean_s, cond_var_s^0.5, log=T)
            #        mh_ratio = mh_ratio - dnorm(x_c_mcmc[iobs,id], cond_mean_c, cond_var_c^0.5, log=T)
            #        mh_ratio = mh_ratio - 2*log(r_mcmc[iobs, id])

            #        alpha_mh = min(1, exp(mh_ratio))
            #        alpha_r[iobs, id] = alpha_r[iobs, id] + alpha_mh



            #        if(runif(1,0,1)<alpha_mh)
            #        {
            #            r_mcmc[iobs, id] = r_prop
            #            x_s_mcmc[iobs,id] = x_s_prop[iobs,id]
            #            x_c_mcmc[iobs,id] = x_c_prop[iobs,id]
            #        }else{
            #            x_s_prop[iobs,id] = x_s_mcmc[iobs,id]
            #            x_c_prop[iobs,id] = x_c_mcmc[iobs,id]
            #        }
            #    }
            # }
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
                    cond_mean_c <- kappa_mcmc[id] - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - kappa_mcmc[-id]))
                    cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

                    cond_mean <- matrix(c(cond_mean_c, cond_mean_s), nrow = 2)

                    uvec <- matrix(c(cos(theta[iobs, id] - mu_mcmc[id]), sin(theta[iobs, id] - mu_mcmc[id])), nrow = 2)


                    A <- t(uvec) %*% (inv_cond_sigma) %*% uvec
                    B <- t(uvec) %*% (inv_cond_sigma) %*% cond_mean


                    v1sim <- runif(1, 0, exp(-0.5 * A * (r_mcmc[iobs, id] - B / A)^2))
                    v2sim <- runif(1, 0, 1)
                    rho1 <- B / A + max(-B / A, -sqrt(-2 * log(v1sim) / A))
                    rho2 <- B / A + sqrt(-2 * log(v1sim) / A)

                    r_mcmc[iobs, id] <- ((rho2^2 - rho1^2) * v2sim + rho1^2)^(1 / 2)

                    x_s_mcmc[iobs, id] <- r_mcmc[iobs, id] * sin(theta[iobs, id] - mu_mcmc[id])
                    x_c_mcmc[iobs, id] <- r_mcmc[iobs, id] * cos(theta[iobs, id] - mu_mcmc[id])
                }
            }

            ### sigma


            nu <- par_sigma_adapt + d + 1
            prop_sigma <- rinvwishart(nu, sigma_s_mcmc * (nu - (d + 1)))

            test_sigma(prop_sigma)

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
                    mh_ratio <- mh_ratio + (-0.5 * c(log_det_c_prop$modulus) - 0.5 * t(x_c_mcmc[iobs, ] - kappa_mcmc) %*% lambda_c_prop %*% (x_c_mcmc[iobs, ] - kappa_mcmc))
                    mh_ratio <- mh_ratio + (-0.5 * c(log_det_s_prop$modulus) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_prop %*% (x_s_mcmc[iobs, ]))

                    mh_ratio <- mh_ratio - (-0.5 * c(log_det_c_mcmc$modulus) - 0.5 * t(x_c_mcmc[iobs, ] - kappa_mcmc) %*% lambda_c_mcmc %*% (x_c_mcmc[iobs, ] - kappa_mcmc))
                    mh_ratio <- mh_ratio - (-0.5 * c(log_det_s_mcmc$modulus) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_mcmc %*% (x_s_mcmc[iobs, ]))
                }
                # prior
                mh_ratio <- mh_ratio + dinvwishart(sigma_s_prop, prior_sigma_nu, prior_sigma_psi, log = T)
                mh_ratio <- mh_ratio - dinvwishart(sigma_s_mcmc, prior_sigma_nu, prior_sigma_psi, log = T)

                ## proposal
                mh_ratio <- mh_ratio - dinvwishart(sigma_s_prop, nu, sigma_s_mcmc * (nu - (d + 1)), log = T)
                mh_ratio <- mh_ratio + dinvwishart(sigma_s_mcmc, nu, sigma_s_prop * (nu - (d + 1)), log = T)

                alpha_sigma <- alpha_sigma + min(1, exp(mh_ratio))

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
                # alpha_r = alpha_r/adapt_batch
                alpha_sigma <- alpha_sigma / adapt_batch
                # print(alpha_sigma)
                # print(par_sigma_adapt)
                if ((sum_iter < burnin)) {
                    for (id in 1:d)
                    {
                        sd_mu[id] <- exp(log(sd_mu[id]) + adapt_a / (adapt_b + sum_iter) * (alpha_mu[id] - adapt_alpha_target))
                        alpha_mu[id] <- 0
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
        kappa_out[imcmc, ] <- kappa_mcmc
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
            app <- dmvnorm(r_mcmc[iobs, tf_missig_out[iobs, ]] * cos(theta[iobs, tf_missig_out[iobs, ]] - mu_mcmc[tf_missig_out[iobs, ]]), kappa_mcmc[tf_missig_out[iobs, ]], sigma_c_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)
            app <- app + dmvnorm(r_mcmc[iobs, tf_missig_out[iobs, ]] * sin(theta[iobs, tf_missig_out[iobs, ]] - mu_mcmc[tf_missig_out[iobs, ]]), vec_zero[tf_missig_out[iobs, ]], sigma_s_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)

            for (id in 1:d)
            {
                if (tf_missig_out[iobs, id] == T) {
                    app <- app + log(r_mcmc[iobs, id])
                }
            }


            sum_log_dens_data[iobs] <- sum_log_dens_data[iobs] + app
            sum_dens_data[iobs] <- sum_dens_data[iobs] + exp(app)
        }
    }
    waic_llpd <- sum(log(sum_dens_data[iobs] / sample_to_save))

    p_waic <- 2 * sum(log(sum_dens_data[iobs] / sample_to_save) - sum_log_dens_data[iobs] / sample_to_save)
    return(list(mu_out = mu_out, kappa_out = kappa_out, sigma_s_out = sigma_s_out, sigma_c_out = sigma_c_out, r_out = r_out, missig_out = missig_out, waic = 2 * (waic_llpd - p_waic)))
}

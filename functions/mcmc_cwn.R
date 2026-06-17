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
  na_index = list(NA),
  do_only_ESS = TRUE,
  n_test_sigma = 10,
  type_ess = 1,
  do_ind = FALSE
) {
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
    r_mcmc <- r_init
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

    chol_sigma_c_mcmc <- cholesky(sigma_c_mcmc)
    chol_sigma_s_mcmc <- cholesky(sigma_s_mcmc)

    log_det_c_mcmc <- 2 * sum(log(diag(chol_sigma_c_mcmc)))
    log_det_s_mcmc <- 2 * sum(log(diag(chol_sigma_s_mcmc)))

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
    if (there_are_na) {
        for (id in 1:d) {
            if (length(na_index[[id]]) > 0) {
                tf_missig_out[na_index[[id]], id] <- FALSE
            }
        }
    }

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

    sum_sq_log_dens_data <- rep(0, n)

    log_sum_dens_data <- rep(-Inf, n)

    n_samp_cop <- 10000
    samp_wc <- matrix(NA, nrow = n_samp_cop, ncol = d)
    ess_acc_sigma <- rep(NA, burn_thin + thin * (sample_to_save - 1))
    counts_ess <- rep(NA, burn_thin + thin * (sample_to_save - 1))
    metropolis_acc_sigma <- rep(NA, burn_thin + thin * (sample_to_save - 1))
    # EES
    n_par_sigma <- d * (d + 1) / 2
    psi_sigma_inv <- solve(prior_sigma_psi)
    chol_psi_sigma_inv <- t(cholesky(psi_sigma_inv))
    inv_chol_psi_sigma_inv <- solve(chol_psi_sigma_inv)


    chol_lambda_s_mcmc <- t(cholesky(lambda_s_mcmc))
    idx_lower <- lower.tri(chol_lambda_s_mcmc, diag = FALSE)
    L_prop_raw <- matrix(0, nrow = d, ncol = d)
    f_r_ess <- matrix(0, nrow = n, ncol = d)
    f_prime_r_ess <- matrix(0, nrow = n, ncol = d)
    v_r_ess <- matrix(0, nrow = n, ncol = d)
    L_prop_raw <- matrix(0, nrow = d, ncol = d)
    molt_ident <- diag(d)
    chol_lambda_s_prop <- matrix(0, nrow = d, ncol = d)


    ### bartlett
    # par_nu_post <- prior_sigma_nu + n
    # par_xi <- par_nu_post + 1 - (1:d)
    # par_psi_post_s <- prior_sigma_psi
    # par_psi_post_s <- par_psi_post_s + crossprod(x_s_mcmc)

    # chol_par_psi_post_s <- t(cholesky(chol2inv(cholesky(par_psi_post_s))))
    # inv_chol_par_psi_post_s <- solve(chol_par_psi_post_s)


    # chol_lambda_s_mcmc <- t(cholesky(lambda_s_mcmc))
    # init


    # Xc_cent <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")
    # Xs <- x_s_mcmc
    # Sc <- crossprod(x_c_mcmc)
    # mu_fc <- log(par_xi) / 2
    # var_fc <- (2 / (par_xi * 4))
    # sd_fc <- sqrt(var_fc)
    ## Ss <- crossprod(x_s_mcmc)

    ## bartlett_mcmc <- inv_chol_par_psi_post_s %*% chol_lambda_s_mcmc
    # bartlett_mcmc <- forwardsolve(
    #    chol_par_psi_post_s,
    #    chol_lambda_s_mcmc,
    #    upper.tri = FALSE
    # )
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

                        # mh_ratio <- mh_ratio + log(func_d_wc(theta_prop, mu_mcmc[id], rho_mcmc[id]))
                        # mh_ratio <- mh_ratio - log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]))

                        mh_ratio <- mh_ratio + func_logd_wc(theta_prop, mu_mcmc[id], rho_mcmc[id])
                        mh_ratio <- mh_ratio - func_logd_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id])

                        if (is.na(mh_ratio)) {
                            print("mh_ratio is NA missing")
                            mh_ratio <- log(0)
                        }

                        if (log(runif(1, 0, 1)) < (mh_ratio)) {
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

            #### mu
            # NOTE: OLD
            # x_s_prop <- x_s_mcmc
            # x_c_prop <- x_c_mcmc
            # theta_cop_prop <- theta_cop
            # for (id in 1:d)
            # {
            #    cond_var_c <- 1 / lambda_c_mcmc[id, id]
            #    cond_var_s <- 1 / lambda_s_mcmc[id, id]

            #    mu_prop <- rnorm(1, mu_mcmc[id], sd_mu[id])


            #    for (iobs in 1:n)
            #    {
            #        theta_cop_prop[iobs, id] <- 2 * pi * func_cdf_wc(theta[iobs, id] - mu_prop, 0, rho_mcmc[id])
            #    }
            #    x_s_prop[, id] <- r_mcmc[, id] * sin(theta_cop_prop[, id])
            #    x_c_prop[, id] <- r_mcmc[, id] * cos(theta_cop_prop[, id])

            #    mh_ratio <- 0
            #    for (iobs in 1:n)
            #    {
            #        cond_mean_c <- 0 - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - 0))
            #        cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

            #        mh_ratio <- mh_ratio + dnorm(x_c_prop[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
            #        mh_ratio <- mh_ratio + dnorm(x_s_prop[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

            #        mh_ratio <- mh_ratio - dnorm(x_c_mcmc[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
            #        mh_ratio <- mh_ratio - dnorm(x_s_mcmc[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

            #        mh_ratio <- mh_ratio + log(func_d_wc(theta[iobs, id], mu_prop, rho_mcmc[id]))

            #        mh_ratio <- mh_ratio - log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]))
            #    }

            #    mh_ratio <- mh_ratio + dnorm(mu_prop, prior_mu_mean[id], prior_mu_var[id]^0.5, log = T)
            #    mh_ratio <- mh_ratio - dnorm(mu_mcmc[id], prior_mu_mean[id], prior_mu_var[id]^0.5, log = T)

            #    alpha_mh <- min(1, exp(mh_ratio))


            #    # if(id == 1)
            #    # {
            #    #    print(alpha_mh)
            #    #    print(sd_mu[id])
            #    # }
            #    if (is.na(alpha_mh)) {
            #        print("mh_ratio is NA mu")
            #        alpha_mh <- 0
            #    }
            #    alpha_mu[id] <- alpha_mu[id] + alpha_mh
            #    if (runif(1, 0, 1) < alpha_mh) {
            #        mu_mcmc[id] <- mu_prop
            #        x_s_mcmc[, id] <- x_s_prop[, id]
            #        x_c_mcmc[, id] <- x_c_prop[, id]
            #        theta_cop[, id] <- theta_cop_prop[, id]
            #    } else {
            #        x_s_prop[, id] <- x_s_mcmc[, id]
            #        x_c_prop[, id] <- x_c_mcmc[, id]
            #        theta_cop_prop[, id] <- theta_cop[, id]
            #    }
            # }
            # NOTE new
            x_s_prop <- x_s_mcmc
            x_c_prop <- x_c_mcmc
            theta_cop_prop <- theta_cop

            for (id in 1:d)
            {
                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]
                sd_c <- sqrt(cond_var_c)
                sd_s <- sqrt(cond_var_s)

                mu_prop <- rnorm(1, mu_mcmc[id], sd_mu[id])

                theta_cop_prop[, id] <- 2 * pi * func_cdf_wc(theta[, id] - mu_prop, 0, rho_mcmc[id])

                x_s_prop[, id] <- r_mcmc[, id] * sin(theta_cop_prop[, id])
                x_c_prop[, id] <- r_mcmc[, id] * cos(theta_cop_prop[, id])

                cond_mean_c_vec <- -cond_var_c * as.vector(
                    x_c_mcmc[, -id, drop = FALSE] %*% lambda_c_mcmc[id, -id]
                )

                cond_mean_s_vec <- -cond_var_s * as.vector(
                    x_s_mcmc[, -id, drop = FALSE] %*% lambda_s_mcmc[id, -id]
                )

                mh_ratio <- sum(dnorm(x_c_prop[, id], cond_mean_c_vec, sd_c, log = TRUE)) +
                    sum(dnorm(x_s_prop[, id], cond_mean_s_vec, sd_s, log = TRUE)) -
                    sum(dnorm(x_c_mcmc[, id], cond_mean_c_vec, sd_c, log = TRUE)) -
                    sum(dnorm(x_s_mcmc[, id], cond_mean_s_vec, sd_s, log = TRUE))

                mh_ratio <- mh_ratio +
                    sum(func_logd_wc(theta[, id], mu_prop, rho_mcmc[id])) -
                    sum(func_logd_wc(theta[, id], mu_mcmc[id], rho_mcmc[id]))

                mh_ratio <- mh_ratio +
                    dnorm(mu_prop, prior_mu_mean[id], sqrt(prior_mu_var[id]), log = TRUE) -
                    dnorm(mu_mcmc[id], prior_mu_mean[id], sqrt(prior_mu_var[id]), log = TRUE)

                alpha_mh <- min(1, exp(mh_ratio))

                if (is.na(alpha_mh)) {
                    print("mh_ratio is NA mu")
                    alpha_mh <- 0
                }

                alpha_mu[id] <- alpha_mu[id] + alpha_mh

                if (log(runif(1)) < mh_ratio) {
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
            # NOTE old
            # x_s_prop <- x_s_mcmc
            # x_c_prop <- x_c_mcmc
            # theta_cop_prop <- theta_cop
            # for (id in 1:d)
            # {
            #    cond_var_c <- 1 / lambda_c_mcmc[id, id]
            #    cond_var_s <- 1 / lambda_s_mcmc[id, id]

            #    rho_mcmc_app <- log(rho_mcmc[id] / (1 - rho_mcmc[id]))
            #    rho_prop_app <- rnorm(1, rho_mcmc_app, sd_rho[id])
            #    # rho_prop_app = rho_mcmc_app
            #    rho_prop <- exp(rho_prop_app) / (1 + exp(rho_prop_app))


            #    for (iobs in 1:n)
            #    {
            #        theta_cop_prop[iobs, id] <- 2 * pi * func_cdf_wc(theta[iobs, id] - mu_mcmc[id], 0, rho_prop)
            #    }
            #    x_s_prop[, id] <- r_mcmc[, id] * sin(theta_cop_prop[, id])
            #    x_c_prop[, id] <- r_mcmc[, id] * cos(theta_cop_prop[, id])

            #    mh_ratio <- 0
            #    for (iobs in 1:n)
            #    {
            #        cond_mean_c <- 0 - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - 0))
            #        cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

            #        mh_ratio <- mh_ratio + dnorm(x_c_prop[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
            #        mh_ratio <- mh_ratio + dnorm(x_s_prop[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

            #        mh_ratio <- mh_ratio - dnorm(x_c_mcmc[iobs, id], cond_mean_c, cond_var_c^0.5, log = T)
            #        mh_ratio <- mh_ratio - dnorm(x_s_mcmc[iobs, id], cond_mean_s, cond_var_s^0.5, log = T)

            #        mh_ratio <- mh_ratio + log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_prop))

            #        mh_ratio <- mh_ratio - log(func_d_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]))
            #    }

            #    mh_ratio <- mh_ratio + (dbeta(rho_prop, prior_rho_a[id], prior_rho_b[id], log = T) + rho_prop_app - 2 * log(1 + exp(rho_prop_app)))
            #    mh_ratio <- mh_ratio - (dbeta(rho_mcmc[id], prior_rho_a[id], prior_rho_b[id], log = T) + rho_mcmc_app - 2 * log(1 + exp(rho_mcmc_app)))
            #    # print(round(mh_ratio,4))

            #    alpha_mh <- min(1, exp(mh_ratio))

            #    # rho = exp(x)/(1+exp(x))

            #    # f_x(x) = f_r(r) dr/dx
            #    if (is.na(alpha_mh)) {
            #        print("mh_ratio is NA rho")
            #        alpha_mh <- 0
            #    }
            #    alpha_rho[id] <- alpha_rho[id] + alpha_mh
            #    if (runif(1, 0, 1) < alpha_mh) {
            #        rho_mcmc[id] <- rho_prop
            #        x_s_mcmc[, id] <- x_s_prop[, id]
            #        x_c_mcmc[, id] <- x_c_prop[, id]
            #        theta_cop[, id] <- theta_cop_prop[, id]
            #    } else {
            #        x_s_prop[, id] <- x_s_mcmc[, id]
            #        x_c_prop[, id] <- x_c_mcmc[, id]
            #        theta_cop_prop[, id] <- theta_cop[, id]
            #    }
            # }
            # NOTE new
            x_s_prop <- x_s_mcmc
            x_c_prop <- x_c_mcmc
            theta_cop_prop <- theta_cop

            for (id in 1:d)
            {
                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]
                sd_c <- sqrt(cond_var_c)
                sd_s <- sqrt(cond_var_s)

                rho_mcmc_app <- log(rho_mcmc[id] / (1 - rho_mcmc[id]))
                rho_prop_app <- rnorm(1, rho_mcmc_app, sd_rho[id])
                rho_prop <- plogis(rho_prop_app) # exp(rho_prop_app) / (1 + exp(rho_prop_app))

                theta_cop_prop[, id] <- 2 * pi * func_cdf_wc(theta[, id] - mu_mcmc[id], 0, rho_prop)

                x_s_prop[, id] <- r_mcmc[, id] * sin(theta_cop_prop[, id])
                x_c_prop[, id] <- r_mcmc[, id] * cos(theta_cop_prop[, id])

                cond_mean_c_vec <- -cond_var_c * as.vector(
                    x_c_mcmc[, -id, drop = FALSE] %*% lambda_c_mcmc[id, -id]
                )

                cond_mean_s_vec <- -cond_var_s * as.vector(
                    x_s_mcmc[, -id, drop = FALSE] %*% lambda_s_mcmc[id, -id]
                )

                mh_ratio <- sum(dnorm(x_c_prop[, id], cond_mean_c_vec, sd_c, log = TRUE)) +
                    sum(dnorm(x_s_prop[, id], cond_mean_s_vec, sd_s, log = TRUE)) -
                    sum(dnorm(x_c_mcmc[, id], cond_mean_c_vec, sd_c, log = TRUE)) -
                    sum(dnorm(x_s_mcmc[, id], cond_mean_s_vec, sd_s, log = TRUE))

                mh_ratio <- mh_ratio +
                    sum(func_logd_wc(theta[, id], mu_mcmc[id], rho_prop)) -
                    sum(func_logd_wc(theta[, id], mu_mcmc[id], rho_mcmc[id]))

                mh_ratio <- mh_ratio +
                    (dbeta(rho_prop, prior_rho_a[id], prior_rho_b[id], log = TRUE) +
                        rho_prop_app - 2 * log(1 + exp(rho_prop_app))) -
                    (dbeta(rho_mcmc[id], prior_rho_a[id], prior_rho_b[id], log = TRUE) +
                        rho_mcmc_app - 2 * log(1 + exp(rho_mcmc_app)))

                alpha_mh <- min(1, exp(mh_ratio))

                if (is.na(alpha_mh)) {
                    print("mh_ratio is NA rho")
                    alpha_mh <- 0
                }

                alpha_rho[id] <- alpha_rho[id] + alpha_mh

                if (log(runif(1)) < mh_ratio) {
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
            # NOTE old
            # for (id in 1:d)
            # {
            #    x_s_prop <- x_s_mcmc
            #    x_c_prop <- x_c_mcmc

            #    cond_var_c <- 1 / lambda_c_mcmc[id, id]
            #    cond_var_s <- 1 / lambda_s_mcmc[id, id]
            #    cond_sigma <- diag(c(1 / lambda_c_mcmc[id, id], 1 / lambda_s_mcmc[id, id]))
            #    inv_cond_sigma <- solve(cond_sigma)
            #    for (iobs in 1:n)
            #    {
            #        cond_mean_c <- 0 - cond_var_c * sum(lambda_c_mcmc[id, -id] * (x_c_mcmc[iobs, -id] - 0))
            #        cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id] * x_s_mcmc[iobs, -id])

            #        cond_mean <- matrix(c(cond_mean_c, cond_mean_s), nrow = 2)

            #        uvec <- matrix(c(cos(theta_cop[iobs, id] - 0), sin(theta_cop[iobs, id] - 0)), nrow = 2)


            #        A <- t(uvec) %*% (inv_cond_sigma) %*% uvec
            #        B <- t(uvec) %*% (inv_cond_sigma) %*% cond_mean


            #        v1sim <- runif(1, 0, exp(-0.5 * A * (r_mcmc[iobs, id] - B / A)^2))
            #        v2sim <- runif(1, 0, 1)
            #        rho1 <- B / A + max(-B / A, -sqrt(-2 * log(v1sim) / A))
            #        rho2 <- B / A + sqrt(-2 * log(v1sim) / A)

            #        r_mcmc[iobs, id] <- ((rho2^2 - rho1^2) * v2sim + rho1^2)^(1 / 2)

            #        x_s_mcmc[iobs, id] <- r_mcmc[iobs, id] * sin(theta_cop[iobs, id] - 0)
            #        x_c_mcmc[iobs, id] <- r_mcmc[iobs, id] * cos(theta_cop[iobs, id] - 0)
            #    }
            # }
            # NOTE new

            for (id in 1:d)
            {
                cond_var_c <- 1 / lambda_c_mcmc[id, id]
                cond_var_s <- 1 / lambda_s_mcmc[id, id]

                prec_c <- lambda_c_mcmc[id, id]
                prec_s <- lambda_s_mcmc[id, id]

                cc <- cos(theta_cop[, id])
                ss <- sin(theta_cop[, id])

                cond_mean_c_vec <- -cond_var_c * as.vector(
                    x_c_mcmc[, -id, drop = FALSE] %*% lambda_c_mcmc[id, -id]
                )

                cond_mean_s_vec <- -cond_var_s * as.vector(
                    x_s_mcmc[, -id, drop = FALSE] %*% lambda_s_mcmc[id, -id]
                )

                A_vec <- prec_c * cc^2 + prec_s * ss^2
                B_vec <- prec_c * cc * cond_mean_c_vec + prec_s * ss * cond_mean_s_vec

                BA_vec <- B_vec / A_vec
                quad_vec <- -0.5 * A_vec * (r_mcmc[, id] - BA_vec)^2

                u1 <- runif(n)
                u2 <- runif(n)

                # v1sim <- u1 * exp(quad_vec)
                # rad <- sqrt(-2 * log(v1sim) / A_vec)

                # rho1 <- BA_vec + pmax(-BA_vec, -rad)
                # rho2 <- BA_vec + rad

                log_v1sim <- log(u1) + quad_vec
                rad <- sqrt(-2 * log_v1sim / A_vec)
                rho1 <- BA_vec + pmax(-BA_vec, -rad)
                rho2 <- BA_vec + rad

                r_mcmc[, id] <- sqrt((rho2^2 - rho1^2) * u2 + rho1^2)

                x_s_mcmc[, id] <- r_mcmc[, id] * ss
                x_c_mcmc[, id] <- r_mcmc[, id] * cc
            }
            ### sigma
            # print(cov(x_s_mcmc[, 1], x_s_mcmc[, 3]))
            # print(c(x_s_mcmc[1, 1:6]))
            # print(theta_cop[1, ])
            # print(r_mcmc[1, 1:6])


            # error("")
            # nu <- par_sigma_adapt + d + 1
            ## prop_sigma <- rinvwishart(nu, sigma_s_mcmc * (nu - (d + 1)))
            ## prop_sigma <- rInvWishart(1, nu, sigma_s_mcmc * (nu - (d + 1)))[, , 1]
            # par_nu_post <- nu + n
            # par_psi_post <- prior_sigma_psi
            # for (iobs in 1:n)
            # {
            #    par_psi_post <- par_psi_post + t(x_s_mcmc[iobs, , drop = F]) %*% (x_s_mcmc[iobs, , drop = F])
            # }
            # ccc <- runif(1, 0, 1)
            # prop_sigma <- ccc * rInvWishart(1, par_nu_post, par_psi_post)[, , 1] + (1 - ccc) * sigma_s_mcmc


            # prop_sigma <- (prop_sigma + t(prop_sigma)) / 2
            ## test_sigma(prop_sigma)

            ## if ((test_sigma(prop_sigma) == TRUE) & (test_sigma2(prop_sigma) == TRUE)) {
            # ess_acc_sigma[sum_iter] <- 0
            # Emat <- matrix(rnorm(d^2, 0, par_sigma_adapt), nrow = d, ncol = d)
            # Emat <- (Emat + t(Emat)) / sqrt(2)
            # diag(Emat) <- rnorm(d, 0, par_sigma_adapt)

            # prop_sigma <- sigma_s_mcmc + Emat
            if (do_ind == FALSE) {
                # SECTION: ESS su FULL CONDITIONAL X_s

                if (do_only_ESS == TRUE) {
                    ess_acc_sigma[sum_iter] <- 0
                    counts_ess[sum_iter] <- 0
                    if (type_ess == 1) {
                        


                        par_nu_post <- prior_sigma_nu + n
                        par_xi <- par_nu_post + 1 - (1:d)
                        par_psi_post_s <- prior_sigma_psi
                        par_psi_post_s <- par_psi_post_s + crossprod(x_s_mcmc)

                        chol_par_psi_post_s <- t(cholesky(chol2inv(cholesky(par_psi_post_s))))
                        # inv_chol_par_psi_post_s <- solve(chol_par_psi_post_s)


                        # chol_lambda_s_mcmc <- t(cholesky(lambda_s_mcmc))
                        # init


                        # Xc_cent <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")
                        # Xs <- x_s_mcmc
                        Sc <- crossprod(x_c_mcmc)
                        mu_fc <- log(par_xi) / 2
                        var_fc <- (2 / (par_xi * 4))
                        sd_fc <- sqrt(var_fc)
                        # Ss <- crossprod(x_s_mcmc)

                        # bartlett_mcmc <- inv_chol_par_psi_post_s %*% chol_lambda_s_mcmc
                        bartlett_mcmc <- forwardsolve(
                            chol_par_psi_post_s,
                            chol_lambda_s_mcmc,
                            upper.tri = FALSE
                        )

                        for (itest in 1:n_test_sigma)
                        {
                            f_m_ess <- bartlett_mcmc[idx_lower]
                            f_c_ess <- log(diag(bartlett_mcmc))

                            # NOTE ESS

                            v_m_ess <- rnorm(n_par_sigma - d, 0, 1)
                            v_c_ess <- rnorm(d, mu_fc, sd_fc)

                            log_diag_mcmc <- f_c_ess
                            u_ess <- runif(1, 0, 1)
                            log_y <- log(u_ess)
                            # log_y <- log_y + (-0.5 * n * log_det_c_mcmc - 0.5 * sum(lambda_c_mcmc * Sc) - 0.5 * n * log_det_s_mcmc - 0.5 * sum(lambda_s_mcmc * Ss))
                            # log_y <- log_y + (-sum(dnorm(log_diag_mcmc, log = T)) + sum(dchisq(exp(2 * log_diag_mcmc), df = prior_sigma_nu + 1 - (1:d), log = TRUE) + log(2) + 2 * log_diag_mcmc))
                            log_y <- log_y + (-0.5 * n * log_det_c_mcmc - 0.5 * sum(lambda_c_mcmc * Sc))
                            log_y <- log_y + (-sum(dnorm(log_diag_mcmc, mu_fc, sd_fc, log = T)) + sum(dchisq(exp(2 * log_diag_mcmc), df = par_xi, log = TRUE) + log(2) + 2 * log_diag_mcmc))

                            #
                            theta_ess <- runif(1, 0, 2 * pi)
                            theta_min_ess <- theta_ess - 2 * pi
                            theta_max_ess <- theta_ess

                            #
                            max_ess_steps <- 500
                            ess_step <- 0
                            acc <- FALSE

                            while ((acc == FALSE) && (ess_step < max_ess_steps)) {
                                if (!is.finite(theta_min_ess) || !is.finite(theta_max_ess) ||

                                    (theta_max_ess - theta_min_ess) < 1e-12) {
                                    break
                                }
                                cos_t <- cos(theta_ess)
                                sin_t <- sin(theta_ess)
                                ess_step <- ess_step + 1
                                f_prime_m_ess <- f_m_ess * cos_t + v_m_ess * sin_t
                                f_prime_c_ess <- mu_fc + (f_c_ess - mu_fc) * cos_t + (v_c_ess - mu_fc) * sin_t

                                if (min(f_prime_c_ess) < -80 || any(!is.finite(f_prime_c_ess))) {
                                    if (theta_ess < 0) {
                                        theta_min_ess <- theta_ess
                                    } else {
                                        theta_max_ess <- theta_ess
                                    }

                                    if (theta_max_ess <= theta_min_ess) break

                                    theta_ess <- runif(1, theta_min_ess, theta_max_ess)

                                    next
                                }
                                L_prop_raw[] <- 0
                                L_prop_raw[idx_lower] <- f_prime_m_ess
                                diag(L_prop_raw) <- exp(f_prime_c_ess)
                                # chol_lambda_s_prop[idx_lower] <- f_prime_m_ess
                                # diag(chol_lambda_s_prop) <- exp(f_prime_c_ess)
                                chol_lambda_s_prop <- chol_par_psi_post_s %*% L_prop_raw

                                prop_sigma <- tryCatch(
                                    chol2inv(t(chol_lambda_s_prop)),
                                    error = function(e) NULL
                                )
                                if (is.null(prop_sigma) || any(!is.finite(prop_sigma))) {
                                    if (theta_ess < 0) {
                                        theta_min_ess <- theta_ess
                                    } else {
                                        theta_max_ess <- theta_ess
                                    }

                                    if (theta_max_ess <= theta_min_ess) break

                                    theta_ess <- runif(1, theta_min_ess, theta_max_ess)

                                    next
                                }
                                # prop_sigma <- chol2inv(t(chol_lambda_s_prop)) #solve(chol_lambda_s_prop %*% t(chol_lambda_s_prop))

                                prop_sigma <- (prop_sigma + t(prop_sigma)) / 2
                                res_test <- test_sigma_mcmc(prop_sigma)
                                counts_ess[sum_iter] <- counts_ess[sum_iter] + 1
                                if (res_test$ind == TRUE) {
                                    ess_acc_sigma[sum_iter] <- ess_acc_sigma[sum_iter] + 1
                                    sigma_s_prop <- prop_sigma
                                    sigma_c_prop <- abs(prop_sigma)

                                    sigma_s_prop <- (sigma_s_prop + t(sigma_s_prop)) / 2
                                    sigma_c_prop <- (sigma_c_prop + t(sigma_c_prop)) / 2

                                    chol_sigma_c_prop <- res_test$chol_sigma_c
                                    chol_sigma_s_prop <- res_test$chol_sigma_s
                                    lambda_c_prop <- chol2inv(chol_sigma_c_prop)
                                    lambda_s_prop <- chol_lambda_s_prop %*% t(chol_lambda_s_prop)

                                    log_det_c_prop <- 2 * sum(log(diag(chol_sigma_c_prop)))
                                    log_det_s_prop <- 2 * sum(log(diag(chol_sigma_s_prop)))

                                    log_diag_prop <- f_prime_c_ess
                                    # log_d_prop <- (-0.5 * n * log_det_c_prop - 0.5 * sum(lambda_c_prop * Sc) - 0.5 * n * log_det_s_prop - 0.5 * sum(lambda_s_prop * Ss))
                                    # log_d_prop <- log_d_prop + (-sum(dnorm(log_diag_prop, log = T)) + sum(dchisq(exp(2 * log_diag_prop), df = prior_sigma_nu + 1 - (1:d), log = TRUE) + log(2) + 2 * log_diag_prop))
                                    log_d_prop <- (-0.5 * n * log_det_c_prop - 0.5 * sum(lambda_c_prop * Sc))
                                    log_d_prop <- log_d_prop + (-sum(dnorm(log_diag_prop, mu_fc, sd_fc, log = T)) + sum(dchisq(exp(2 * log_diag_prop), df = par_xi, log = TRUE) + log(2) + 2 * log_diag_prop))

                                    if (is.finite(log_d_prop) && log_d_prop > log_y) {
                                        acc <- TRUE


                                        sigma_s_mcmc <- sigma_s_prop
                                        sigma_c_mcmc <- sigma_c_prop

                                        lambda_s_mcmc <- lambda_s_prop
                                        lambda_c_mcmc <- lambda_c_prop

                                        log_det_c_mcmc <- log_det_c_prop
                                        log_det_s_mcmc <- log_det_s_prop

                                        chol_lambda_s_mcmc <- chol_lambda_s_prop

                                        bartlett_mcmc <- L_prop_raw
                                    }
                                } else {
                                    
                                }
                                if (theta_ess < 0) {
                                    theta_min_ess <- theta_ess
                                } else {
                                    theta_max_ess <- theta_ess
                                }
                                theta_ess <- runif(1, theta_min_ess, theta_max_ess)
                            }
                            if (!acc) {}
                        }
                    }
                    # SECTION: ESS su PRIOR
                    if (type_ess == 2) {
                        error("Type 2 ESS not implemented")   
                    }
                    # SECTION: ESS su PRIOR con R
                    if (type_ess == 3) {
                        error("Type 3 ESS not implemented")
                    }
                    ess_acc_sigma[sum_iter] <- ess_acc_sigma[sum_iter] / counts_ess[sum_iter]
                } else {
                    metropolis_acc_sigma[sum_iter] <- 0
                    for (itest in 1:n_test_sigma)
                    {
                        nu_adapt <- par_sigma_adapt + d +1
                        #prop_sigma <- rWishart(1, nu_adapt, sigma_s_mcmc / nu_adapt)[, , 1]
                        prop_sigma <- rInvWishart(1, nu_adapt, (nu_adapt- d - 1 ) *sigma_s_mcmc  )[, , 1]

                        prop_sigma <- (prop_sigma + t(prop_sigma)) / 2

                        # print("Test Sigma")
                        res_test <- test_sigma_mcmc(prop_sigma)
                        if (res_test$ind == TRUE) {
                            # print("A")
                            metropolis_acc_sigma[sum_iter] <- metropolis_acc_sigma[sum_iter] + 1
                            sigma_s_prop <- prop_sigma
                            sigma_c_prop <- abs(prop_sigma)

                            sigma_s_prop <- (sigma_s_prop + t(sigma_s_prop)) / 2
                            sigma_c_prop <- (sigma_c_prop + t(sigma_c_prop)) / 2

                            chol_sigma_c_prop <- res_test$chol_sigma_c
                            chol_sigma_s_prop <- res_test$chol_sigma_s
                            lambda_c_prop <- chol2inv(chol_sigma_c_prop)
                            lambda_s_prop <- chol2inv(chol_sigma_s_prop)

                            # log_det_c_mcmc <-
                            log_det_c_prop <- 2 * sum(log(diag(chol_sigma_c_prop)))

                            # log_det_s_mcmc <-
                            log_det_s_prop <- 2 * sum(log(diag(chol_sigma_s_prop)))


                            mh_ratio <- 0
                            # for (iobs in 1:n)
                            # {
                            #    mh_ratio <- mh_ratio + (-0.5 * c(log_det_c_prop) - 0.5 * t(x_c_mcmc[iobs, ] - 0) %*% lambda_c_prop %*% (x_c_mcmc[iobs, ] - 0))
                            #    mh_ratio <- mh_ratio + (-0.5 * c(log_det_s_prop) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_prop %*% (x_s_mcmc[iobs, ]))

                            #    mh_ratio <- mh_ratio - (-0.5 * c(log_det_c_mcmc) - 0.5 * t(x_c_mcmc[iobs, ] - 0) %*% lambda_c_mcmc %*% (x_c_mcmc[iobs, ] - 0))
                            #    mh_ratio <- mh_ratio - (-0.5 * c(log_det_s_mcmc) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_mcmc %*% (x_s_mcmc[iobs, ]))
                            # }
                            Sc <- crossprod(x_c_mcmc)
                            Ss <- crossprod(x_s_mcmc)

                            mh_ratio <- mh_ratio + (
                                -0.5 * n * log_det_c_prop - 0.5 * sum(lambda_c_prop * Sc) -
                                    0.5 * n * log_det_s_prop - 0.5 * sum(lambda_s_prop * Ss) +
                                    0.5 * n * log_det_c_mcmc + 0.5 * sum(lambda_c_mcmc * Sc) +
                                    0.5 * n * log_det_s_mcmc + 0.5 * sum(lambda_s_mcmc * Ss)
                            )

                            # prior
                            mh_ratio <- mh_ratio + dInvWishart(sigma_s_prop, prior_sigma_nu, prior_sigma_psi, log = T)
                            mh_ratio <- mh_ratio - dInvWishart(sigma_s_mcmc, prior_sigma_nu, prior_sigma_psi, log = T)

                            ### proposal
                            mh_ratio <- mh_ratio - (dInvWishart(prop_sigma, nu_adapt, (nu_adapt- d - 1 )  *sigma_s_mcmc , log = T))
                            mh_ratio <- mh_ratio + (dInvWishart(sigma_s_mcmc, nu_adapt, (nu_adapt- d - 1 ) *prop_sigma , log = T))


                            if (is.na(mh_ratio)) {
                                print("New NA alpha")
                                mh_ratio <- -Inf
                            }
                            alpha_sigma <- alpha_sigma + min(1, exp(mh_ratio)) / n_test_sigma
                            if (log(runif(1, 0, 1)) < mh_ratio) {
                                # print("ACC")
                                sigma_s_mcmc <- sigma_s_prop
                                sigma_c_mcmc <- sigma_c_prop

                                lambda_s_mcmc <- lambda_s_prop
                                lambda_c_mcmc <- lambda_c_prop

                                log_det_c_mcmc <- log_det_c_prop
                                log_det_s_mcmc <- log_det_s_prop

                                chol_lambda_s_mcmc <- t(cholesky(lambda_s_mcmc))
                            }
                        }
                    }
                    metropolis_acc_sigma[sum_iter] <- metropolis_acc_sigma[sum_iter] / n_test_sigma
                }
            } else {
                sigma_s_mcmc <- diag(1, d)
                sigma_c_mcmc <- diag(1, d)
                lambda_s_mcmc <- diag(1, d)
                lambda_c_mcmc <- diag(1, d)
                log_det_c_mcmc <- 0
                log_det_s_mcmc <- 0
            }


            # if (do_only_ESS == TRUE) {} else {

            # }
            ## ESS SIGMA AND RHO
            # if (do_only_ESS == TRUE) {
            #    # PRIOR/SULL CONDITIONAL
            #    par_nu_post <- prior_sigma_nu + n
            #    par_xi <- par_nu_post + 1 - (1:d)
            #    par_psi_post_s <- prior_sigma_psi
            #    for (iobs in 1:n)
            #    {
            #        par_psi_post_s <- par_psi_post_s + t(x_s_mcmc[iobs, , drop = F]) %*% (x_s_mcmc[iobs, , drop = F])
            #    }
            #    chol_par_psi_post_s <- t(cholesky(chol2inv(cholesky(par_psi_post_s))))
            #    inv_chol_par_psi_post_s <- solve(chol_par_psi_post_s)


            #    chol_lambda_s_mcmc <- t(cholesky(lambda_s_mcmc))

            #    # init
            #    chol_lambda_s_prop <- matrix(0, nrow = d, ncol = d)


            #    # Xc_cent <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")
            #    # Xs <- x_s_mcmc
            #    Sc <- crossprod(x_c_mcmc)
            #    # Ss <- crossprod(Xs)

            #    bartlett_mcmc <- inv_chol_par_psi_post_s %*% chol_lambda_s_mcmc

            #    # NOTE:initial variables
            #    f_m_ess <- bartlett_mcmc[idx_lower]
            #    f_c_ess <- log(diag(bartlett_mcmc))
            #    f_r_ess <- log(r_mcmc)

            #    # SIMULATION FORM "PRIOR"
            #    # simul R
            #    sss_r <- 1
            #    mu_r_ess <- rep(0.0, d)
            #    var_r_ess <- rep(sss_r, d)
            #    for (id in 1:d)
            #    {
            #        v_r_ess[, id] <- rnorm(n, mu_r_ess[id], var_r_ess[id]^0.5)
            #    }
            #    for (id in 1:d)
            #    {
            #        x_c_prop[, id] <- exp(v_r_ess[, id]) * cos(theta[, id])
            #        x_s_prop[, id] <- exp(v_r_ess[, id]) * sin(theta[, id])
            #    }

            #    mu_fc <- log(par_xi) / 2
            #    var_fc <- (2 / (par_xi * 4))
            #    v_m_ess <- rnorm(n_par_sigma - d, 0, 1)
            #    v_c_ess <- rnorm(d, mu_fc, (var_fc)^0.5)
            #    log_diag_mcmc <- f_c_ess

            #    # NOTE: Likelihood mcmc
            #    u_ess <- runif(1, 0, 1)
            #    log_y <- log(u_ess)
            #    log_y <- log_y + (-0.5 * n * log_det_c_mcmc - 0.5 * sum(lambda_c_mcmc * Sc))
            #    log_y <- log_y + (-sum(dnorm(log_diag_mcmc, mu_fc, (var_fc)^0.5, log = T)) + sum(dchisq(exp(2 * log_diag_mcmc), df = par_xi, log = TRUE) + log(2) + 2 * log_diag_mcmc))
            #    for (id in 1:d)
            #    {
            #        log_y <- log_y - sum(dnorm(f_r_ess[, id], mu_r_ess[id], var_r_ess[id]^0.5, log = T))
            #    }
            #    log_y <- log_y + sum(f_r_ess)


            #    # log_y <- log_y + (-0.5 * n * log_det_c_mcmc - 0.5 * sum(lambda_c_mcmc * Sc) - 0.5 * n * log_det_s_mcmc - 0.5 * sum(lambda_s_mcmc * Ss))
            #    # log_y <- log_y + (-sum(dnorm(log_diag_mcmc, log = T)) + sum(dchisq(exp(2 * log_diag_mcmc), df = prior_sigma_nu + 1 - (1:d), log = TRUE) + log(2) + 2 * log_diag_mcmc))

            #    #
            #    theta_ess <- runif(1, 0, 2 * pi)
            #    theta_min_ess <- theta_ess - 2 * pi
            #    theta_max_ess <- theta_ess

            #    #
            #    acc <- FALSE
            #    ess_acc_sigma[sum_iter] <- 0
            #    while (acc == FALSE) {
            #        f_prime_m_ess <- f_m_ess * cos(theta_ess) + v_m_ess * sin(theta_ess)
            #        f_prime_c_ess <- mu_fc + (f_c_ess - mu_fc) * cos(theta_ess) + (v_c_ess - mu_fc) * sin(theta_ess)
            #        for (id in 1:d)
            #        {
            #            f_prime_r_ess[, id] <- mu_r_ess[id] + (f_r_ess[, id] - mu_r_ess[id]) * cos(theta_ess) + (v_r_ess[, id] - mu_r_ess[id]) * sin(theta_ess)
            #        }

            #        if (min(f_prime_c_ess) < -45) {
            #            # reject this angle

            #            if (theta_ess < 0) {
            #                theta_min_ess <- theta_ess
            #            } else {
            #                theta_max_ess <- theta_ess
            #            }

            #            theta_ess <- runif(1, theta_min_ess, theta_max_ess)

            #            next
            #        }
            #        L_prop_raw[idx_lower] <- f_prime_m_ess
            #        diag(L_prop_raw) <- exp(f_prime_c_ess)
            #        # chol_lambda_s_prop[idx_lower] <- f_prime_m_ess
            #        # diag(chol_lambda_s_prop) <- exp(f_prime_c_ess)
            #        chol_lambda_s_prop <- chol_par_psi_post_s %*% L_prop_raw

            #        prop_sigma <- tryCatch(
            #            chol2inv(t(chol_lambda_s_prop)),
            #            error = function(e) NULL
            #        )
            #        if (is.null(prop_sigma)) {
            #            if (theta_ess < 0) {
            #                theta_min_ess <- theta_ess
            #            } else {
            #                theta_max_ess <- theta_ess
            #            }

            #            theta_ess <- runif(1, theta_min_ess, theta_max_ess)

            #            next
            #        }

            #        prop_sigma <- (prop_sigma + t(prop_sigma)) / 2
            #        res_test <- test_sigma_mcmc(prop_sigma)
            #        if (res_test$ind == TRUE) {
            #            sigma_s_prop <- prop_sigma
            #            sigma_c_prop <- abs(prop_sigma)

            #            sigma_s_prop <- (sigma_s_prop + t(sigma_s_prop)) / 2
            #            sigma_c_prop <- (sigma_c_prop + t(sigma_c_prop)) / 2

            #            chol_sigma_c_prop <- res_test$chol_sigma_c
            #            chol_sigma_s_prop <- res_test$chol_sigma_s
            #            lambda_c_prop <- chol2inv(chol_sigma_c_prop)
            #            lambda_s_prop <- chol_lambda_s_prop %*% t(chol_lambda_s_prop)

            #            log_det_c_prop <- 2 * sum(log(diag(chol_sigma_c_prop)))
            #            log_det_s_prop <- 2 * sum(log(diag(chol_sigma_s_prop)))

            #            log_diag_prop <- f_prime_c_ess
            #            for (id in 1:d)
            #            {
            #                x_c_prop[, id] <- exp(f_prime_r_ess[, id]) * cos(theta[, id])
            #                x_s_prop[, id] <- exp(f_prime_r_ess[, id]) * sin(theta[, id])
            #            }

            #            # Xc_cent <- sweep(x_c_prop, 2, kappa_mcmc, FUN = "-")
            #            # Xs <- x_s_mcmc
            #            Sc_prop <- crossprod(x_c_prop)
            #            # log_d_prop <- (-0.5 * n * log_det_c_prop - 0.5 * sum(lambda_c_prop * Sc) - 0.5 * n * log_det_s_prop - 0.5 * sum(lambda_s_prop * Ss))
            #            # log_d_prop <- log_d_prop + (-sum(dnorm(log_diag_prop, log = T)) + sum(dchisq(exp(2 * log_diag_prop), df = prior_sigma_nu + 1 - (1:d), log = TRUE) + log(2) + 2 * log_diag_prop))
            #            log_d_prop <- (-0.5 * n * log_det_c_prop - 0.5 * sum(lambda_c_prop * Sc_prop))
            #            log_d_prop <- log_d_prop + (-sum(dnorm(log_diag_prop, mu_fc, (var_fc)^0.5, log = T)) + sum(dchisq(exp(2 * log_diag_prop), df = par_xi, log = TRUE) + log(2) + 2 * log_diag_prop))
            #            for (id in 1:d)
            #            {
            #                log_d_prop <- log_d_prop - sum(dnorm(f_prime_r_ess[, id], mu_r_ess[id], var_r_ess[id]^0.5, log = T))
            #            }
            #            log_d_prop <- log_d_prop + sum(f_prime_r_ess)

            #            if (is.finite(log_d_prop) && log_d_prop > log_y) {
            #                acc <- TRUE


            #                sigma_s_mcmc <- sigma_s_prop
            #                sigma_c_mcmc <- sigma_c_prop

            #                lambda_s_mcmc <- lambda_s_prop
            #                lambda_c_mcmc <- lambda_c_prop

            #                log_det_c_mcmc <- log_det_c_prop
            #                log_det_s_mcmc <- log_det_s_prop

            #                chol_lambda_s_mcmc <- chol_lambda_s_prop

            #                r_mcmc <- exp(f_prime_r_ess)
            #                x_c_mcmc <- x_c_prop
            #                x_s_mcmc <- x_s_prop
            #            }
            #        } else {
            #            ess_acc_sigma[sum_iter] <- ess_acc_sigma[sum_iter] + 1
            #        }
            #        if (theta_ess < 0) {
            #            theta_min_ess <- theta_ess
            #        } else {
            #            theta_max_ess <- theta_ess
            #        }
            #        theta_ess <- runif(1, theta_min_ess, theta_max_ess)
            #    }


            #    # lambda_chol_p <- matrix(0, nrow = d, ncol = d)
            #    # L[lower.tri(L, diag = TRUE)] <- v
            # }


            ### update of the adaptive parameters
            if ((sum_iter %% adapt_batch == 0)) {
                alpha_mu <- alpha_mu / adapt_batch
                alpha_rho <- alpha_rho / adapt_batch
                # alpha_r = alpha_r/adapt_batch
                alpha_sigma <- alpha_sigma / adapt_batch
                # print(cbind(alpha_rho,sd_rho))
                # print(par_sigma_adapt)
                if ((sum_iter < (burnin * 100.0))) {
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
                    if ((do_only_ESS == FALSE) & (do_ind == FALSE)) {
                        #print(c(par_sigma_adapt, alpha_sigma))
                        par_sigma_adapt <- exp(log(par_sigma_adapt) - adapt_a / (adapt_b + sum_iter) * (alpha_sigma - adapt_alpha_target))
                    }

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
        # for (iobs in 1:n)
        # {
        #    for (id in 1:d)
        #    {
        #        theta_zeta[id] <- 2 * pi * func_cdf_wc(theta[iobs, id] - mu_mcmc[id], 0, rho_mcmc[id])
        #        x_s_zeta[id] <- r_mcmc[iobs, id] * sin(theta_zeta[id] - 0)
        #        x_c_zeta[id] <- r_mcmc[iobs, id] * cos(theta_zeta[id] - 0)
        #    }
        #    app <- dmvnorm(x_c_zeta[tf_missig_out[iobs, ]], vec_zero[tf_missig_out[iobs, ]], sigma_c_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)

        #    app <- app + dmvnorm(x_s_zeta[tf_missig_out[iobs, ]], vec_zero[tf_missig_out[iobs, ]], sigma_s_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)

        #    for (id in 1:d)
        #    {
        #        if (tf_missig_out[iobs, id] == T) {
        #            app <- app + func_logd_wc(theta[iobs, id], mu_mcmc[id], rho_mcmc[id]) + log(r_mcmc[iobs, id]) + log(2 * pi)
        #        }
        #    }


        #    sum_log_dens_data[iobs] <- sum_log_dens_data[iobs] + app
        #    sum_sq_log_dens_data[iobs] <- sum_sq_log_dens_data[iobs] + app^2
        #    log_sum_dens_data[iobs] <- logsumexp2(log_sum_dens_data[iobs], app)
        # }
    }


    # waic_llpd <- sum(log_sum_dens_data - log(sample_to_save))

    # mean_log_dens <- sum_log_dens_data / sample_to_save
    # var_log_dens <- sum_sq_log_dens_data / sample_to_save - mean_log_dens^2

    # p_waic <- sum(var_log_dens)

    return(list(mu_out = mu_out, rho_out = rho_out, sigma_s_out = sigma_s_out, sigma_c_out = sigma_c_out, r_out = r_out, missig_out = missig_out, waic = NA, ess_acc_sigma = list(ess_acc_sigma, metropolis_acc_sigma, counts_ess)))
}

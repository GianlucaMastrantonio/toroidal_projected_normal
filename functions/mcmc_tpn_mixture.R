log_tpn_latent_density_i <- function(
  iobs, k, theta, r_mcmc,
  mu_mcmc, kappa_mcmc,
  lambda_c_mcmc, lambda_s_mcmc,
  log_det_c_mcmc, log_det_s_mcmc,
  d
) {
    xc <- r_mcmc[iobs, ] * cos(theta[iobs, ] - mu_mcmc[1, , k])
    xs <- r_mcmc[iobs, ] * sin(theta[iobs, ] - mu_mcmc[1, , k])

    xc_cent <- xc - kappa_mcmc[1, , k]

    out <- -d * log(2 * pi)

    out <- out -
        0.5 * log_det_c_mcmc[k] -
        0.5 * as.numeric(t(xc_cent) %*% lambda_c_mcmc[, , k] %*% xc_cent)

    out <- out -
        0.5 * log_det_s_mcmc[k] -
        0.5 * as.numeric(t(xs) %*% lambda_s_mcmc[, , k] %*% xs)

    out
}

update_z_mfm_tpn <- function(
  z_mcmc,
  nvec,
  kappa_mcmc,
  lambda_c_mcmc,
  lambda_s_mcmc,
  log_det_c_mcmc,
  log_det_s_mcmc,
  logV,
  gamma_mfm,
  Kmax,
  d,
  theta,
  r_mcmc,
  mu_mcmc
) {
    n <- length(z_mcmc)

    for (iobs in 1:n) {
        old_k <- z_mcmc[iobs]
        nvec[old_k] <- nvec[old_k] - 1

        occupied <- which(nvec > 0)
        empty <- which(nvec == 0)
        t_minus_i <- length(occupied)

        logw <- rep(-Inf, Kmax)

        for (k in occupied) {
            logw[k] <- log(nvec[k] + gamma_mfm) +
                log_tpn_latent_density_i(
                    iobs = iobs,
                    k = k,
                    theta = theta,
                    r_mcmc = r_mcmc,
                    mu_mcmc = mu_mcmc,
                    kappa_mcmc = kappa_mcmc,
                    lambda_c_mcmc = lambda_c_mcmc,
                    lambda_s_mcmc = lambda_s_mcmc,
                    log_det_c_mcmc = log_det_c_mcmc,
                    log_det_s_mcmc = log_det_s_mcmc,
                    d = d
                )
        }

        if (length(empty) > 0 && t_minus_i < Kmax) {
            k_new <- empty[1]

            if (t_minus_i == 0) {
                log_ratio_V <- 0
            } else {
                log_ratio_V <- logV[t_minus_i + 1] - logV[t_minus_i]
            }

            logw[k_new] <- log(gamma_mfm) +
                log_ratio_V +
                log_tpn_latent_density_i(
                    iobs = iobs,
                    k = k_new,
                    theta = theta,
                    r_mcmc = r_mcmc,
                    mu_mcmc = mu_mcmc,
                    kappa_mcmc = kappa_mcmc,
                    lambda_c_mcmc = lambda_c_mcmc,
                    lambda_s_mcmc = lambda_s_mcmc,
                    log_det_c_mcmc = log_det_c_mcmc,
                    log_det_s_mcmc = log_det_s_mcmc,
                    d = d
                )
        }

        prob <- exp(logw - logsumexp(logw))
        new_k <- sample(seq_len(Kmax), size = 1, prob = prob)

        z_mcmc[iobs] <- new_k
        nvec[new_k] <- nvec[new_k] + 1
    }

    list(z_mcmc = z_mcmc, nvec = nvec)
}
update_z_finite_tpn <- function(
  z_mcmc,
  nvec,
  kappa_mcmc,
  lambda_c_mcmc,
  lambda_s_mcmc,
  log_det_c_mcmc,
  log_det_s_mcmc,
  gamma,
  K,
  d,
  theta,
  r_mcmc,
  mu_mcmc
) {
    n <- length(z_mcmc)

    for (iobs in 1:n) {
        old_k <- z_mcmc[iobs]
        nvec[old_k] <- nvec[old_k] - 1

        logw <- rep(-Inf, K)

        for (k in 1:K) {
            logw[k] <- log(nvec[k] + gamma) +
                log_tpn_latent_density_i(
                    iobs = iobs,
                    k = k,
                    theta = theta,
                    r_mcmc = r_mcmc,
                    mu_mcmc = mu_mcmc,
                    kappa_mcmc = kappa_mcmc,
                    lambda_c_mcmc = lambda_c_mcmc,
                    lambda_s_mcmc = lambda_s_mcmc,
                    log_det_c_mcmc = log_det_c_mcmc,
                    log_det_s_mcmc = log_det_s_mcmc,
                    d = d
                )
        }

        prob <- exp(logw - logsumexp(logw))
        new_k <- sample(seq_len(K), size = 1, prob = prob)

        z_mcmc[iobs] <- new_k
        nvec[new_k] <- nvec[new_k] + 1
    }

    list(z_mcmc = z_mcmc, nvec = nvec)
}


mcmc_tpn_mixture <- function(
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
  na_index = list(NA),
  do_only_ESS = TRUE,
  n_test_sigma = 10,
  type_ess = 1,
  do_ind = FALSE,
  Kmax = 10,
  gamma_mfm = 1
) {
    log_pK <- dpois(1:Kmax, lambda = 1, log = TRUE)
    log_pK <- log_pK - logsumexp(log_pK)

    d <- dim(theta)[2]
    n <- dim(theta)[1]
    sample_to_save <- round((iterations - burnin) / thin)

    # This containt the posterior sample that we are going to save
    mu_out <- array(NA, c(sample_to_save, d, Kmax))
    kappa_out <- array(NA, c(sample_to_save, d, Kmax))
    sigma_s_out <- array(NA, c(sample_to_save, d^2, Kmax))
    sigma_c_out <- array(NA, c(sample_to_save, d^2, Kmax))
    r_out <- array(NA, c(sample_to_save, n, d))

    # these objects contains the corrent value of the parameters
    mu_mcmc <- array(NA, c(1, d, Kmax))
    kappa_mcmc <- array(NA, c(1, d, Kmax))
    sigma_s_mcmc <- array(NA, c(d, d, Kmax))
    sigma_c_mcmc <- array(NA, c(d, d, Kmax))
    lambda_c_mcmc <- array(NA, c(d, d, Kmax))
    lambda_s_mcmc <- array(NA, c(d, d, Kmax))
    r_mcmc <- matrix(abs(rnorm(n * d)), nrow = n, ncol = d)
    # for(iobs in 1:n)
    # {
    #    for(id in 1:d)
    #    {
    #        r_mcmc[iobs, id] = (x_c[iobs, id]^2+x_s[iobs, id]^2)^0.5
    #    }
    # }
    z_mcmc <- sample(1:min(5, Kmax), n, replace = TRUE)
    nvec <- rep(0, Kmax)
    for (k in 1:Kmax)
    {
        nvec[k] <- sum(z_mcmc == k)
    }
    logV <- precompute_logV(
        n = n,
        H = Kmax,
        gamma = gamma_mfm,
        log_pK = log_pK
    )
    z_out <- array(NA, c(sample_to_save, n))
    # the linear variables
    x_s_mcmc <- matrix(1, nrow = n, ncol = d)
    x_c_mcmc <- matrix(1, nrow = n, ncol = d)

    for (k in 1:Kmax)
    {
        mu_mcmc[1, , k] <- mu_init
        kappa_mcmc[1, , k] <- kappa_init
    }
    r_mcmc[] <- r_init

    for (id in 1:d)
    {
        x_s_mcmc[, id] <- r_mcmc[, id] * sin(theta[, id] - mu_mcmc[1, id, z_mcmc])
        x_c_mcmc[, id] <- r_mcmc[, id] * cos(theta[, id] - mu_mcmc[1, id, z_mcmc])
    }
    x_s_prop <- x_s_mcmc
    x_c_prop <- x_c_mcmc

    # the adaptive part of the Metropolis
    # sd_r <- matrix(1, nrow = n, ncol = d)
    # alpha_r <- matrix(0, nrow = n, ncol = d)
    sd_mu <- matrix(1, nrow = d, ncol = Kmax) * sd_mu_scal
    alpha_mu <- matrix(0, nrow = d, ncol = Kmax)
    # alpha_sigma <- 0


    # other objects that containt the current value
    for (k in 1:Kmax)
    {
        sigma_s_mcmc[, , k] <- (sigma_init + t(sigma_init)) / 2 # i did this because soimethimes, the matrices are not exactly simmetrical
        sigma_c_mcmc[, , k] <- abs(sigma_s_mcmc[, , k])
    }

    chol_sigma_c_mcmc <- chol_sigma_s_mcmc <- sigma_s_mcmc
    log_det_c_mcmc <- log_det_s_mcmc <- rep(NA, Kmax)
    lambda_c_mcmc <- lambda_s_mcmc <- sigma_s_mcmc
    for (k in 1:Kmax)
    {
        chol_sigma_c_mcmc[, , k] <- cholesky(sigma_c_mcmc[, , k])
        chol_sigma_s_mcmc[, , k] <- cholesky(sigma_s_mcmc[, , k])

        log_det_c_mcmc[k] <- 2 * sum(log(diag(chol_sigma_c_mcmc[, , k])))
        log_det_s_mcmc[k] <- 2 * sum(log(diag(chol_sigma_s_mcmc[, , k])))


        lambda_c_mcmc[, , k] <- chol2inv(cholesky(sigma_c_mcmc[, , k]))
        lambda_s_mcmc[, , k] <- chol2inv(cholesky(sigma_s_mcmc[, , k]))
    }


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
    # vec_zero <- matrix(0, nrow = d)
    # sum_log_dens_data <- rep(0, n)
    # sum_dens_data <- rep(0, n)
    # sum_log_dens_data <- rep(0, n)
    # sum_sq_log_dens_data <- rep(0, n)
    # log_sum_dens_data <- rep(-Inf, n)
    ess_acc_sigma <- rep(NA, burn_thin + thin * (sample_to_save - 1))
    counts_ess <- rep(NA, burn_thin + thin * (sample_to_save - 1))
    metropolis_acc_sigma <- rep(NA, burn_thin + thin * (sample_to_save - 1))
    # prop_prec_sigma_mcmc <- 0.5
    # sigma_iw_mcmc <- sigma_s_mcmc
    # sigma_iw_prop <- sigma_s_mcmc
    # prop_prec_sigma_prop <- 0.5
    # EES
    n_par_sigma <- d * (d + 1) / 2
    psi_sigma_inv <- solve(prior_sigma_psi)
    # chol_psi_sigma_inv <- t(cholesky(psi_sigma_inv))
    # inv_chol_psi_sigma_inv <- solve(chol_psi_sigma_inv)


    # chol_lambda_s_mcmc <- t(cholesky(lambda_s_mcmc))
    idx_lower <- lower.tri(matrix(0, d, d), diag = FALSE)
    # f_r_ess <- matrix(0, nrow = n, ncol = d)
    # f_prime_r_ess <- matrix(0, nrow = n, ncol = d)
    # v_r_ess <- matrix(0, nrow = n, ncol = d)
    L_prop_raw <- matrix(0, nrow = d, ncol = d)

    chol_lambda_s_prop <- matrix(0, nrow = d, ncol = d)


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

            # for (k in which(nvec == 0)) {
            #    mu_mcmc[1, , k] <- rnorm(d, prior_mu_mean, sqrt(prior_mu_var))

            #    kappa_mcmc[1, , k] <- rtruncnorm(
            #        d,
            #        a = 0,
            #        b = Inf,
            #        mean = prior_kappa_mean,
            #        sd = sqrt(prior_kappa_var)
            #    )
            # }
            # res_z <- update_z_finite_tpn(
            #    z_mcmc = z_mcmc,
            #    nvec = nvec,
            #    kappa_mcmc = kappa_mcmc,
            #    lambda_c_mcmc = lambda_c_mcmc,
            #    lambda_s_mcmc = lambda_s_mcmc,
            #    log_det_c_mcmc = log_det_c_mcmc,
            #    log_det_s_mcmc = log_det_s_mcmc,
            #    gamma = gamma_mfm,
            #    K = Kmax,
            #    d = d,
            #    theta = theta,
            #    r_mcmc = r_mcmc,
            #    mu_mcmc = mu_mcmc
            # )

            res_z <- update_z_mfm_tpn(
                z_mcmc = z_mcmc,
                nvec = nvec,
                kappa_mcmc = kappa_mcmc,
                lambda_c_mcmc = lambda_c_mcmc,
                lambda_s_mcmc = lambda_s_mcmc,
                log_det_c_mcmc = log_det_c_mcmc,
                log_det_s_mcmc = log_det_s_mcmc,
                logV = logV,
                gamma_mfm = gamma_mfm,
                Kmax = Kmax,
                d = d,
                theta = theta,
                r_mcmc = r_mcmc,
                mu_mcmc = mu_mcmc
            )
            z_mcmc <- res_z$z_mcmc
            nvec <- res_z$nvec
            for (id in 1:d) {
                x_s_mcmc[, id] <- r_mcmc[, id] *
                    sin(theta[, id] - mu_mcmc[1, id, z_mcmc])

                x_c_mcmc[, id] <- r_mcmc[, id] *
                    cos(theta[, id] - mu_mcmc[1, id, z_mcmc])
            }
            #### missing
            if (there_are_na) {
                for (id in 1:d)
                {
                    for (iobs in na_index[[id]])
                    {
                        k <- z_mcmc[iobs]
                        cond_var_c <- 1 / lambda_c_mcmc[id, id, k]
                        cond_var_s <- 1 / lambda_s_mcmc[id, id, k]

                        cond_mean_c <- kappa_mcmc[1, id, k] - cond_var_c * sum(lambda_c_mcmc[id, -id, k] * (x_c_mcmc[iobs, -id] - kappa_mcmc[1, -id, k]))
                        cond_mean_s <- -cond_var_s * sum(lambda_s_mcmc[id, -id, k] * x_s_mcmc[iobs, -id])

                        x_c_mcmc[iobs, id] <- rnorm(1, cond_mean_c, cond_var_c^0.5)
                        x_s_mcmc[iobs, id] <- rnorm(1, cond_mean_s, cond_var_s^0.5)

                        r_mcmc[iobs, id] <- sqrt((x_c_mcmc[iobs, id]^2 + x_s_mcmc[iobs, id]^2))
                        theta[iobs, id] <- atan2(x_s_mcmc[iobs, id], x_c_mcmc[iobs, id]) + mu_mcmc[1, id, k]
                    }
                }
            }


            # NOTE NEW updata
            # x_s_prop <- x_s_mcmc
            # x_c_prop <- x_c_mcmc
            # xc_centered <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")
            # for (id in 1:d)
            # {
            #    cond_var_c <- 1 / lambda_c_mcmc[id, id]
            #    cond_var_s <- 1 / lambda_s_mcmc[id, id]
            #    sd_c <- sqrt(cond_var_c)
            #    sd_s <- sqrt(cond_var_s)

            #    mu_prop <- rnorm(1, mu_mcmc[id], sd_mu[id])

            #    x_s_prop[, id] <- r_mcmc[, id] * sin(theta[, id] - mu_prop)
            #    x_c_prop[, id] <- r_mcmc[, id] * cos(theta[, id] - mu_prop)


            #    mh_ratio <- 0

            #    # xc_centered <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")

            #    cond_mean_c_vec <- kappa_mcmc[id] - cond_var_c * (xc_centered[, -id, drop = FALSE] %*% lambda_c_mcmc[id, -id])
            #    cond_mean_s_vec <- -cond_var_s * (x_s_mcmc[, -id, drop = FALSE] %*% lambda_s_mcmc[id, -id])

            #    mh_ratio <- sum(dnorm(x_c_prop[, id], cond_mean_c_vec, sd_c, log = TRUE)) +
            #        sum(dnorm(x_s_prop[, id], cond_mean_s_vec, sd_s, log = TRUE)) -
            #        sum(dnorm(x_c_mcmc[, id], cond_mean_c_vec, sd_c, log = TRUE)) -
            #        sum(dnorm(x_s_mcmc[, id], cond_mean_s_vec, sd_s, log = TRUE))

            #    mh_ratio <- mh_ratio + dnorm(mu_prop, prior_mu_mean[id], prior_mu_var[id]^0.5, log = T)
            #    mh_ratio <- mh_ratio - dnorm(mu_mcmc[id], prior_mu_mean[id], prior_mu_var[id]^0.5, log = T)

            #    alpha_mh <- min(1, exp(mh_ratio))
            #    alpha_mu[id] <- alpha_mu[id] + alpha_mh

            #    if (log(runif(1, 0, 1)) < mh_ratio) {
            #        mu_mcmc[id] <- mu_prop
            #        x_s_mcmc[, id] <- x_s_prop[, id]
            #        x_c_mcmc[, id] <- x_c_prop[, id]
            #        xc_centered[, id] <- x_c_mcmc[, id] - kappa_mcmc[id]
            #    } else {
            #        x_s_prop[, id] <- x_s_mcmc[, id]
            #        x_c_prop[, id] <- x_c_mcmc[, id]
            #    }
            # }


            x_s_prop <- x_s_mcmc
            x_c_prop <- x_c_mcmc

            for (k in 1:Kmax) {
                ind_k <- which(z_mcmc == k)
                nk <- length(ind_k)

                if (nk > 0) {
                    xc_centered_k <- sweep(
                        x_c_mcmc[ind_k, , drop = FALSE],
                        2,
                        kappa_mcmc[1, , k],
                        FUN = "-"
                    )
                }

                for (id in 1:d) {
                    mu_prop <- rnorm(1, mu_mcmc[1, id, k], sd_mu[id, k])

                    mh_ratio <-
                        dnorm(mu_prop, prior_mu_mean[id], sqrt(prior_mu_var[id]), log = TRUE) -
                        dnorm(mu_mcmc[1, id, k], prior_mu_mean[id], sqrt(prior_mu_var[id]), log = TRUE)

                    if (nk > 0) {
                        cond_var_c <- 1 / lambda_c_mcmc[id, id, k]
                        cond_var_s <- 1 / lambda_s_mcmc[id, id, k]

                        sd_c <- sqrt(cond_var_c)
                        sd_s <- sqrt(cond_var_s)

                        x_s_prop[ind_k, id] <- r_mcmc[ind_k, id] *
                            sin(theta[ind_k, id] - mu_prop)

                        x_c_prop[ind_k, id] <- r_mcmc[ind_k, id] *
                            cos(theta[ind_k, id] - mu_prop)

                        cond_mean_c_vec <- kappa_mcmc[1, id, k] -
                            cond_var_c * as.vector(
                                xc_centered_k[, -id, drop = FALSE] %*%
                                    lambda_c_mcmc[id, -id, k]
                            )

                        cond_mean_s_vec <- -cond_var_s * as.vector(
                            x_s_mcmc[ind_k, -id, drop = FALSE] %*%
                                lambda_s_mcmc[id, -id, k]
                        )

                        mh_ratio <- mh_ratio +
                            sum(dnorm(x_c_prop[ind_k, id], cond_mean_c_vec, sd_c, log = TRUE)) +
                            sum(dnorm(x_s_prop[ind_k, id], cond_mean_s_vec, sd_s, log = TRUE)) -
                            sum(dnorm(x_c_mcmc[ind_k, id], cond_mean_c_vec, sd_c, log = TRUE)) -
                            sum(dnorm(x_s_mcmc[ind_k, id], cond_mean_s_vec, sd_s, log = TRUE))
                    }

                    alpha_mh <- min(1, exp(mh_ratio))
                    alpha_mu[id, k] <- alpha_mu[id, k] + alpha_mh

                    if (log(runif(1)) < mh_ratio) {
                        mu_mcmc[1, id, k] <- mu_prop

                        if (nk > 0) {
                            x_s_mcmc[ind_k, id] <- x_s_prop[ind_k, id]
                            x_c_mcmc[ind_k, id] <- x_c_prop[ind_k, id]
                            xc_centered_k[, id] <- x_c_mcmc[ind_k, id] - kappa_mcmc[1, id, k]
                        }
                    } else {
                        if (nk > 0) {
                            x_s_prop[ind_k, id] <- x_s_mcmc[ind_k, id]
                            x_c_prop[ind_k, id] <- x_c_mcmc[ind_k, id]
                        }
                    }
                }
            }

            #### kappa

            # NOTE: New
            # xc_centered <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")
            # for (id in 1:d)
            # {
            #    cond_var_c <- 1 / lambda_c_mcmc[id, id]

            #    adj_vec <- -cond_var_c * as.vector(
            #        xc_centered[, -id, drop = FALSE] %*% lambda_c_mcmc[id, -id]
            #    )

            #    sum_app_x <- sum(x_c_mcmc[, id] - adj_vec)

            #    par_2 <- 1 / (n / cond_var_c + 1 / prior_kappa_var[id])
            #    par_1 <- par_2 * (sum_app_x / cond_var_c + prior_kappa_mean[id] / prior_kappa_var[id])

            #    kappa_mcmc[id] <- rtruncnorm(1, a = 0, b = Inf, mean = par_1, sd = sqrt(par_2))
            #    xc_centered[, id] <- x_c_mcmc[, id] - kappa_mcmc[id]
            # }

            for (k in 1:Kmax) {
                ind_k <- which(z_mcmc == k)
                nk <- length(ind_k)

                for (id in 1:d) {
                    cond_var_c <- 1 / lambda_c_mcmc[id, id, k]

                    if (nk > 0) {
                        xc_centered_k <- sweep(
                            x_c_mcmc[ind_k, , drop = FALSE],
                            2,
                            kappa_mcmc[1, , k],
                            FUN = "-"
                        )

                        adj_vec <- -cond_var_c * as.vector(
                            xc_centered_k[, -id, drop = FALSE] %*%
                                lambda_c_mcmc[id, -id, k]
                        )

                        sum_app_x <- sum(x_c_mcmc[ind_k, id] - adj_vec)

                        par_2 <- 1 / (nk / cond_var_c + 1 / prior_kappa_var[id])
                        par_1 <- par_2 * (
                            sum_app_x / cond_var_c +
                                prior_kappa_mean[id] / prior_kappa_var[id]
                        )

                        kappa_mcmc[1, id, k] <- rtruncnorm(
                            1,
                            a = 0,
                            b = Inf,
                            mean = par_1,
                            sd = sqrt(par_2)
                        )
                    } else {
                        kappa_mcmc[1, id, k] <- rtruncnorm(
                            1,
                            a = 0,
                            b = Inf,
                            mean = prior_kappa_mean[id],
                            sd = sqrt(prior_kappa_var[id])
                        )
                    }
                }
            }
            ## mu and k
            # NOTE: è la full condition del vettore delle medie 2d-variato
            for (k in which(nvec > 0)) {
                ind_k <- which(z_mcmc == k)
                nk <- length(ind_k)

                mu_prop <- rep(NA, d)
                kappa_prop <- rep(NA, d)

                x_c_prop_k <- r_mcmc[ind_k, , drop = FALSE] * cos(theta[ind_k, , drop = FALSE])
                x_s_prop_k <- r_mcmc[ind_k, , drop = FALSE] * sin(theta[ind_k, , drop = FALSE])

                var_p_c <- sigma_c_mcmc[, , k] / nk
                var_p_s <- sigma_s_mcmc[, , k] / nk

                mu_p_c <- lambda_c_mcmc[, , k] %*% matrix(colSums(x_c_prop_k), ncol = 1)
                mu_p_s <- lambda_s_mcmc[, , k] %*% matrix(colSums(x_s_prop_k), ncol = 1)

                mu_p_c <- var_p_c %*% mu_p_c
                mu_p_s <- var_p_s %*% mu_p_s

                prop_s <- mu_p_s + t(cholesky(var_p_s)) %*% rnorm(d)
                prop_c <- mu_p_c + t(cholesky(var_p_c)) %*% rnorm(d)

                for (id in 1:d) {
                    mu_prop[id] <- atan2(prop_s[id], prop_c[id])
                    kappa_prop[id] <- sqrt(prop_s[id]^2 + prop_c[id]^2)
                }

                mh_ratio <- 0

                mh_ratio <- mh_ratio +
                    sum(log(kappa_mcmc[1, , k])) -
                    sum(log(kappa_prop))

                mh_ratio <- mh_ratio +
                    sum(dnorm(kappa_prop, prior_kappa_mean, sqrt(prior_kappa_var), log = TRUE)) -
                    sum(dnorm(kappa_mcmc[1, , k], prior_kappa_mean, sqrt(prior_kappa_var), log = TRUE))

                mh_ratio <- mh_ratio +
                    sum(dnorm(mu_prop, prior_mu_mean, sqrt(prior_mu_var), log = TRUE)) -
                    sum(dnorm(mu_mcmc[1, , k], prior_mu_mean, sqrt(prior_mu_var), log = TRUE))

                if (log(runif(1)) < mh_ratio) {
                    mu_mcmc[1, , k] <- mu_prop
                    kappa_mcmc[1, , k] <- kappa_prop

                    x_c_mcmc[ind_k, ] <- r_mcmc[ind_k, , drop = FALSE] *
                        cos(theta[ind_k, , drop = FALSE] -
                            matrix(mu_mcmc[1, , k], nrow = nk, ncol = d, byrow = TRUE))

                    x_s_mcmc[ind_k, ] <- r_mcmc[ind_k, , drop = FALSE] *
                        sin(theta[ind_k, , drop = FALSE] -
                            matrix(mu_mcmc[1, , k], nrow = nk, ncol = d, byrow = TRUE))
                }
            }
            ##### r


            for (k in which(nvec > 0)) {
                ind_k <- which(z_mcmc == k)
                nk <- length(ind_k)

                xc_centered_k <- sweep(
                    x_c_mcmc[ind_k, , drop = FALSE],
                    2,
                    kappa_mcmc[1, , k],
                    FUN = "-"
                )

                for (id in 1:d) {
                    cond_var_c <- 1 / lambda_c_mcmc[id, id, k]
                    cond_var_s <- 1 / lambda_s_mcmc[id, id, k]

                    prec_c <- lambda_c_mcmc[id, id, k]
                    prec_s <- lambda_s_mcmc[id, id, k]

                    ang <- theta[ind_k, id] - mu_mcmc[1, id, k]
                    cc <- cos(ang)
                    ss <- sin(ang)

                    cond_mean_c_vec <- kappa_mcmc[1, id, k] -
                        cond_var_c * as.vector(
                            xc_centered_k[, -id, drop = FALSE] %*%
                                lambda_c_mcmc[id, -id, k]
                        )

                    cond_mean_s_vec <- -cond_var_s * as.vector(
                        x_s_mcmc[ind_k, -id, drop = FALSE] %*%
                            lambda_s_mcmc[id, -id, k]
                    )

                    A_vec <- prec_c * cc^2 + prec_s * ss^2
                    B_vec <- prec_c * cc * cond_mean_c_vec +
                        prec_s * ss * cond_mean_s_vec

                    BA_vec <- B_vec / A_vec
                    quad_vec <- -0.5 * A_vec * (r_mcmc[ind_k, id] - BA_vec)^2

                    u1 <- runif(nk)
                    u2 <- runif(nk)

                    log_v1sim <- log(u1) + quad_vec
                    rad <- sqrt(-2 * log_v1sim / A_vec)

                    rho1 <- BA_vec + pmax(-BA_vec, -rad)
                    rho2 <- BA_vec + rad

                    r_mcmc[ind_k, id] <- sqrt((rho2^2 - rho1^2) * u2 + rho1^2)

                    x_c_mcmc[ind_k, id] <- r_mcmc[ind_k, id] * cc
                    x_s_mcmc[ind_k, id] <- r_mcmc[ind_k, id] * ss

                    xc_centered_k[, id] <- x_c_mcmc[ind_k, id] - kappa_mcmc[1, id, k]
                }
            }

            ### sigma


            # NOTE: new sigma
            # type_ess <- 1
            # print(c(do_ind, do_only_ESS))
            if (do_ind == FALSE) {
                # molt_ident <- diag(d)

                if (do_only_ESS == TRUE) {
                    ess_acc_sigma[sum_iter] <- 0
                    counts_ess[sum_iter] <- 0

                    if (type_ess == 1) {
                        for (k in 1:Kmax) {
                            ind_k <- which(z_mcmc == k)
                            nk <- length(ind_k)

                            par_nu_post <- prior_sigma_nu + nk
                            par_xi <- par_nu_post + 1 - (1:d)

                            if (nk > 0) {
                                Xs_k <- x_s_mcmc[ind_k, , drop = FALSE]
                                Xc_cent_k <- sweep(
                                    x_c_mcmc[ind_k, , drop = FALSE],
                                    2,
                                    kappa_mcmc[1, , k],
                                    FUN = "-"
                                )

                                par_psi_post_s <- prior_sigma_psi + crossprod(Xs_k)
                                Sc <- crossprod(Xc_cent_k)
                            } else {
                                par_psi_post_s <- prior_sigma_psi
                                Sc <- matrix(0, nrow = d, ncol = d)
                            }

                            chol_par_psi_post_s <- t(cholesky(chol2inv(cholesky(par_psi_post_s))))

                            chol_lambda_s_mcmc_k <- t(cholesky(lambda_s_mcmc[, , k]))

                            mu_fc <- log(par_xi) / 2
                            var_fc <- 2 / (par_xi * 4)
                            sd_fc <- sqrt(var_fc)

                            bartlett_mcmc <- forwardsolve(
                                chol_par_psi_post_s,
                                chol_lambda_s_mcmc_k,
                                upper.tri = FALSE
                            )

                            for (itest in 1:n_test_sigma) {
                                f_m_ess <- bartlett_mcmc[idx_lower]
                                f_c_ess <- log(diag(bartlett_mcmc))

                                v_m_ess <- rnorm(n_par_sigma - d, 0, 1)
                                v_c_ess <- rnorm(d, mu_fc, sd_fc)

                                log_diag_mcmc <- f_c_ess

                                log_y <- log(runif(1))
                                log_y <- log_y +
                                    (-0.5 * nk * log_det_c_mcmc[k] -
                                        0.5 * sum(lambda_c_mcmc[, , k] * Sc))

                                log_y <- log_y +
                                    (
                                        -sum(dnorm(log_diag_mcmc, mu_fc, sd_fc, log = TRUE)) +
                                            sum(
                                                dchisq(exp(2 * log_diag_mcmc), df = par_xi, log = TRUE) +
                                                    log(2) +
                                                    2 * log_diag_mcmc
                                            )
                                    )

                                theta_ess <- runif(1, 0, 2 * pi)
                                theta_min_ess <- theta_ess - 2 * pi
                                theta_max_ess <- theta_ess

                                max_ess_steps <- 500
                                ess_step <- 0
                                acc <- FALSE

                                while (!acc && ess_step < max_ess_steps) {
                                    if (
                                        !is.finite(theta_min_ess) ||
                                            !is.finite(theta_max_ess) ||
                                            (theta_max_ess - theta_min_ess) < 1e-12
                                    ) {
                                        break
                                    }

                                    cos_t <- cos(theta_ess)
                                    sin_t <- sin(theta_ess)
                                    ess_step <- ess_step + 1

                                    f_prime_m_ess <- f_m_ess * cos_t + v_m_ess * sin_t
                                    f_prime_c_ess <- mu_fc +
                                        (f_c_ess - mu_fc) * cos_t +
                                        (v_c_ess - mu_fc) * sin_t

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

                                        log_d_prop <- -0.5 * nk * log_det_c_prop -
                                            0.5 * sum(lambda_c_prop * Sc)

                                        log_d_prop <- log_d_prop +
                                            (
                                                -sum(dnorm(log_diag_prop, mu_fc, sd_fc, log = TRUE)) +
                                                    sum(
                                                        dchisq(exp(2 * log_diag_prop), df = par_xi, log = TRUE) +
                                                            log(2) +
                                                            2 * log_diag_prop
                                                    )
                                            )

                                        if (is.finite(log_d_prop) && log_d_prop > log_y) {
                                            acc <- TRUE

                                            sigma_s_mcmc[, , k] <- sigma_s_prop
                                            sigma_c_mcmc[, , k] <- sigma_c_prop

                                            lambda_s_mcmc[, , k] <- lambda_s_prop
                                            lambda_c_mcmc[, , k] <- lambda_c_prop

                                            log_det_c_mcmc[k] <- log_det_c_prop
                                            log_det_s_mcmc[k] <- log_det_s_prop

                                            chol_sigma_c_mcmc[, , k] <- chol_sigma_c_prop
                                            chol_sigma_s_mcmc[, , k] <- chol_sigma_s_prop

                                            bartlett_mcmc <- L_prop_raw
                                        }
                                    }

                                    if (theta_ess < 0) {
                                        theta_min_ess <- theta_ess
                                    } else {
                                        theta_max_ess <- theta_ess
                                    }

                                    theta_ess <- runif(1, theta_min_ess, theta_max_ess)
                                }
                            }
                        }
                    }

                    if (type_ess == 2) {
                        error("Type 2 ESS not implemented")
                    }

                    if (type_ess == 3) {
                        error("Type 3 ESS not implemented")
                    }

                    if (!is.na(counts_ess[sum_iter]) && counts_ess[sum_iter] > 0) {
                        ess_acc_sigma[sum_iter] <- ess_acc_sigma[sum_iter] / counts_ess[sum_iter]
                    }
                } else {
                    #### NO ESS
                    error("NoESS not implemented")
                    # metropolis_acc_sigma[sum_iter] <- 0
                    # for (itest in 1:n_test_sigma)
                    # {
                    #    nu_adapt <- par_sigma_adapt + d + 1
                    #    # prop_sigma <- rWishart(1, nu_adapt, sigma_s_mcmc / nu_adapt)[, , 1]
                    #    prop_sigma <- rInvWishart(1, nu_adapt, (nu_adapt - d - 1) * sigma_s_mcmc)[, , 1]

                    #    prop_sigma <- (prop_sigma + t(prop_sigma)) / 2

                    #    # print("Test Sigma")
                    #    # print(par_sigma_adapt)
                    #    res_test <- test_sigma_mcmc(prop_sigma)
                    #    if (res_test$ind == TRUE) {
                    #        # print("A")
                    #        metropolis_acc_sigma[sum_iter] <- metropolis_acc_sigma[sum_iter] + 1
                    #        sigma_s_prop <- prop_sigma
                    #        sigma_c_prop <- abs(prop_sigma)

                    #        sigma_s_prop <- (sigma_s_prop + t(sigma_s_prop)) / 2
                    #        sigma_c_prop <- (sigma_c_prop + t(sigma_c_prop)) / 2

                    #        chol_sigma_c_prop <- res_test$chol_sigma_c
                    #        chol_sigma_s_prop <- res_test$chol_sigma_s
                    #        lambda_c_prop <- chol2inv(chol_sigma_c_prop)
                    #        lambda_s_prop <- chol2inv(chol_sigma_s_prop)

                    #        # log_det_c_mcmc <-
                    #        log_det_c_prop <- 2 * sum(log(diag(chol_sigma_c_prop)))

                    #        # log_det_s_mcmc <-
                    #        log_det_s_prop <- 2 * sum(log(diag(chol_sigma_s_prop)))


                    #        mh_ratio <- 0
                    #        # for (iobs in 1:n)
                    #        # {
                    #        #    mh_ratio <- mh_ratio + (-0.5 * c(log_det_c_prop) - 0.5 * t(x_c_mcmc[iobs, ] - kappa_mcmc) %*% lambda_c_prop %*% (x_c_mcmc[iobs, ] - kappa_mcmc))
                    #        #    mh_ratio <- mh_ratio + (-0.5 * c(log_det_s_prop) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_prop %*% (x_s_mcmc[iobs, ]))

                    #        #    mh_ratio <- mh_ratio - (-0.5 * c(log_det_c_mcmc) - 0.5 * t(x_c_mcmc[iobs, ] - kappa_mcmc) %*% lambda_c_mcmc %*% (x_c_mcmc[iobs, ] - kappa_mcmc))
                    #        #    mh_ratio <- mh_ratio - (-0.5 * c(log_det_s_mcmc) - 0.5 * t(x_s_mcmc[iobs, ]) %*% lambda_s_mcmc %*% (x_s_mcmc[iobs, ]))
                    #        # }
                    #        Xc_cent <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")
                    #        Xs <- x_s_mcmc
                    #        Sc <- crossprod(Xc_cent)
                    #        Ss <- crossprod(Xs)
                    #        mh_ratio <- mh_ratio + (
                    #            -0.5 * n * log_det_c_prop - 0.5 * sum(lambda_c_prop * Sc) -
                    #                0.5 * n * log_det_s_prop - 0.5 * sum(lambda_s_prop * Ss) +
                    #                0.5 * n * log_det_c_mcmc + 0.5 * sum(lambda_c_mcmc * Sc) +
                    #                0.5 * n * log_det_s_mcmc + 0.5 * sum(lambda_s_mcmc * Ss))


                    #        # print("Data")
                    #        # print(mh_ratio)

                    #        # prior
                    #        # print("Prior")
                    #        mh_ratio <- mh_ratio + dInvWishart(sigma_s_prop, prior_sigma_nu, prior_sigma_psi, log = T)
                    #        mh_ratio <- mh_ratio - dInvWishart(sigma_s_mcmc, prior_sigma_nu, prior_sigma_psi, log = T)
                    #        # print(mh_ratio)
                    #        ### proposal
                    #        # print("Proposal")
                    #        mh_ratio <- mh_ratio - (dInvWishart(prop_sigma, nu_adapt, (nu_adapt - d - 1) * sigma_s_mcmc, log = T))
                    #        mh_ratio <- mh_ratio + (dInvWishart(sigma_s_mcmc, nu_adapt, (nu_adapt - d - 1) * prop_sigma, log = T))
                    #        # print(mh_ratio)

                    #        if (is.na(exp(mh_ratio))) {
                    #            print("New NA alpha")
                    #            mh_ratio <- -Inf
                    #        }
                    #        # print(mh_ratio)
                    #        alpha_sigma <- alpha_sigma + min(1, exp(mh_ratio)) / n_test_sigma
                    #        if (log(runif(1, 0, 1)) < mh_ratio) {
                    #            # print("ACC")
                    #            sigma_s_mcmc <- sigma_s_prop
                    #            sigma_c_mcmc <- sigma_c_prop

                    #            lambda_s_mcmc <- lambda_s_prop
                    #            lambda_c_mcmc <- lambda_c_prop

                    #            log_det_c_mcmc <- log_det_c_prop
                    #            log_det_s_mcmc <- log_det_s_prop

                    #            chol_lambda_s_mcmc <- t(cholesky(lambda_s_mcmc))

                    #            # sigma_iw_mcmc <- sigma_iw_prop
                    #            # prop_prec_sigma_mcmc <- prop_prec_sigma_prop
                    #        }
                    #    }
                    # }
                    # metropolis_acc_sigma[sum_iter] <- metropolis_acc_sigma[sum_iter] / n_test_sigma
                }
            } else {
                for (k in 1:Kmax) {
                    sigma_s_mcmc[, , k] <- diag(1, d)
                    sigma_c_mcmc[, , k] <- diag(1, d)

                    lambda_s_mcmc[, , k] <- diag(1, d)
                    lambda_c_mcmc[, , k] <- diag(1, d)

                    chol_sigma_s_mcmc[, , k] <- diag(1, d)
                    chol_sigma_c_mcmc[, , k] <- diag(1, d)

                    log_det_s_mcmc[k] <- 0
                    log_det_c_mcmc[k] <- 0
                }
            }


            # if (do_only_ESS == TRUE) {
            # lambda_chol_p <- matrix(0, nrow = d, ncol = d)
            # L[lower.tri(L, diag = TRUE)] <- v
            # } else {}
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


            #    Xc_cent <- sweep(x_c_mcmc, 2, kappa_mcmc, FUN = "-")
            #    # Xs <- x_s_mcmc
            #    Sc <- crossprod(Xc_cent)
            #    # Ss <- crossprod(Xs)

            #    bartlett_mcmc <- inv_chol_par_psi_post_s %*% chol_lambda_s_mcmc

            #    # NOTE:initial variables
            #    f_m_ess <- bartlett_mcmc[idx_lower]
            #    f_c_ess <- log(diag(bartlett_mcmc))
            #    f_r_ess <- log(r_mcmc)

            #    # SIMULATION FORM "PRIOR"
            #    # simul R
            #    sss_r <- 10
            #    kappa_safe <- pmax(kappa_mcmc, 1e-7)
            #    kappa2_safe <- kappa_safe^2
            #    mu_r_ess <- log(kappa_safe) - sss_r / (2 * kappa2_safe)
            #    var_r_ess <- sss_r / kappa2_safe
            #    for (id in 1:d)
            #    {
            #        v_r_ess[, id] <- rnorm(n, mu_r_ess[id], var_r_ess[id]^0.5)
            #    }
            #    for (id in 1:d)
            #    {
            #        x_c_prop[, id] <- exp(v_r_ess[, id]) * cos(theta[, id] - mu_mcmc[id])
            #        x_s_prop[, id] <- exp(v_r_ess[, id]) * sin(theta[, id] - mu_mcmc[id])
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
            #                x_c_prop[, id] <- exp(f_prime_r_ess[, id]) * cos(theta[, id] - mu_mcmc[id])
            #                x_s_prop[, id] <- exp(f_prime_r_ess[, id]) * sin(theta[, id] - mu_mcmc[id])
            #            }

            #            Xc_cent <- sweep(x_c_prop, 2, kappa_mcmc, FUN = "-")
            #            # Xs <- x_s_mcmc
            #            Sc_prop <- crossprod(Xc_cent)
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
                # alpha_r = alpha_r/adapt_batch
                # alpha_sigma <- alpha_sigma / adapt_batch
                # print(alpha_sigma)
                # print(par_sigma_adapt)
                if ((sum_iter < (burnin * 100.0))) {
                    for (id in 1:d)
                    {
                        sd_mu[id, ] <- exp(log(sd_mu[id, ]) + adapt_a / (adapt_b + sum_iter) * (alpha_mu[id, ] - adapt_alpha_target))
                        alpha_mu[id, ] <- 0
                        # for(iobs in 1:n)
                        # {
                        #    sd_r[iobs,id] = exp(log(sd_r[iobs,id]) +  adapt_a/(adapt_b+sum_iter)*(alpha_r[iobs,id] - adapt_alpha_target) )
                        #    alpha_r[iobs,id] = 0
                        # }
                    }

                    # sigma
                    # if ((do_only_ESS == FALSE) & (do_ind == FALSE)) {
                    #    # print(c(par_sigma_adapt, alpha_sigma))
                    #    par_sigma_adapt <- exp(log(par_sigma_adapt) - adapt_a / (adapt_b + sum_iter) * (alpha_sigma - adapt_alpha_target))
                    # }

                    # alpha_sigma <- 0
                }
            }
        }
        burn_thin <- thin

        # i save the current values of the chains
        mu_out[imcmc, , ] <- mu_mcmc[1, , ]
        kappa_out[imcmc, , ] <- kappa_mcmc[1, , ]

        for (k in 1:Kmax) {
            sigma_s_out[imcmc, , k] <- c(sigma_s_mcmc[, , k])
            sigma_c_out[imcmc, , k] <- c(sigma_c_mcmc[, , k])
        }

        z_out[imcmc, ] <- z_mcmc
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
        #    app <- dmvnorm(r_mcmc[iobs, tf_missig_out[iobs, ]] * cos(theta[iobs, tf_missig_out[iobs, ]] - mu_mcmc[tf_missig_out[iobs, ]]), kappa_mcmc[tf_missig_out[iobs, ]], sigma_c_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)
        #    app <- app + dmvnorm(r_mcmc[iobs, tf_missig_out[iobs, ]] * sin(theta[iobs, tf_missig_out[iobs, ]] - mu_mcmc[tf_missig_out[iobs, ]]), vec_zero[tf_missig_out[iobs, ]], sigma_s_mcmc[tf_missig_out[iobs, ], tf_missig_out[iobs, ]], log = T)

        #    for (id in 1:d)
        #    {
        #        if (tf_missig_out[iobs, id] == T) {
        #            app <- app + log(r_mcmc[iobs, id])
        #        }
        #    }


        #    # sum_log_dens_data[iobs] <- sum_log_dens_data[iobs] + app
        #    # sum_dens_data[iobs] <- sum_dens_data[iobs] + exp(app)
        #    sum_log_dens_data[iobs] <- sum_log_dens_data[iobs] + app
        #    sum_sq_log_dens_data[iobs] <- sum_sq_log_dens_data[iobs] + app^2
        #    log_sum_dens_data[iobs] <- logsumexp2(log_sum_dens_data[iobs], app)
        # }
    }

    # waic_llpd <- sum(log_sum_dens_data - log(sample_to_save))

    # mean_log_dens <- sum_log_dens_data / sample_to_save
    # var_log_dens <- sum_sq_log_dens_data / sample_to_save - mean_log_dens^2

    # p_waic <- sum(var_log_dens)


    return(list(
        mu_out = mu_out,
        kappa_out = kappa_out,
        sigma_s_out = sigma_s_out,
        sigma_c_out = sigma_c_out,
        r_out = r_out,
        z_out = z_out,
        missig_out = missig_out,
        waic = NA,
        ess_acc_sigma = list(ess_acc_sigma, metropolis_acc_sigma, counts_ess)
    ))
}

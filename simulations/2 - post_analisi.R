# Post-processing script for the simulation study.
# It reads the posterior samples saved by simulations/1 - launch_tpn.R and
# simulations/1 - launch_ctpn.R, groups files by simulation scenario, combines
# chains, and computes convergence and accuracy diagnostics.

# Plot autocorrelation functions from multiple chains on the same panel.
plot_multi_acf <- function(chain, par_name, par_id, max_lag = 100, par_title = NULL) {
  acf_list <- lapply(chain, function(ch) {
    stats::acf(ch[[par_name]][, par_id], lag.max = max_lag, plot = FALSE)
  })

  y <- range(unlist(lapply(acf_list, function(a) a$acf)))

  plot(
    acf_list[[1]]$lag,
    acf_list[[1]]$acf,
    type = "l",
    ylim = y,
    xlab = "Lag",
    ylab = "ACF",
    main <- if (is.null(par_title)) paste(par_name, "par", par_id) else par_title
  )

  for (cc in 2:length(acf_list)) {
    lines(acf_list[[cc]]$lag, acf_list[[cc]]$acf, col = cc)
  }

  abline(h = 0, lty = 2)
}

# Compute univariate ESS for each column and a determinant-based multivariate
# ESS summary. This is used as an additional diagnostic alongside posterior.
multiESS_base2 <- function(x, max_lag = 100, eps = 1e-8) {
  x <- as.matrix(x)

  n <- nrow(x)
  p <- ncol(x)

  Lambda <- stats::cov(x) + diag(eps, p)

  ess <- apply(x, 2, function(z) {
    ac <- as.numeric(stats::acf(
      z,
      lag.max = min(max_lag, length(z) - 1),
      plot = FALSE
    )$acf)

    rho <- ac[-1]
    k <- which(rho < 0)[1]
    if (!is.na(k) && k > 1) rho <- rho[seq_len(k - 1)]
    if (!is.na(k) && k == 1) rho <- numeric(0)

    tau <- max(1, 1 + 2 * sum(rho))
    length(z) / tau
  })

  Sigma <- Lambda
  for (j in seq_len(p)) {
    Sigma[j, j] <- Lambda[j, j] * n / ess[j]
  }

  logdet_Lambda <- as.numeric(determinant(Lambda, logarithm = TRUE)$modulus)
  logdet_Sigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  mess <- n * exp((logdet_Lambda - logdet_Sigma) / p)

  list(
    mESS = mess,
    ESS = ess,
    minESS = min(ess),
    medianESS = median(ess)
  )
}
library(stringr)
library(coda)
library(posterior)

# Select the simulation output files generated with the standard prefix.
name <- " Simulation"
ff <- list.files("simulations/output/")
ff <- ff[startsWith(ff, name)]
ff_cw <- ff[grepl("cwc", ff)]
ff_pn <- ff[grepl("tpn", ff)]

# Parse TPN file names into a data frame of simulation settings.
word_to_remove <- c(name, "do_best_init=", "tpn_simulations_results -  select_d= ", " select_n= ", " select_kappa= ", " select_sigma= ", " select_chain= ", " seed_sigma= ", " seed_data= ", ".Rdata")

pattern <- paste(str_escape(word_to_remove), collapse = "|")
ff_pn_clean <- str_remove_all(ff_pn, pattern)
split_list <- strsplit(ff_pn_clean, " ")
mat_pn <- do.call(rbind, split_list)
data_pn <- as.data.frame(mat_pn)
colnames(data_pn) <- c("do_best_init", "select_d", "select_n", "select_kappa", "select_sigma", "select_chain", "seed_sigma", "seed_data")
i <- 1
data_pn[[i]] <- as.logical(as.character(data_pn[[i]]))
for (i in 2:ncol(data_pn)) {
  data_pn[[i]] <- as.numeric(as.character(data_pn[[i]]))
}
data_pn[, 1] <- FALSE

# Parse copula-model file names into a data frame of simulation settings.
word_to_remove <- c(name, "do_best_init=", "cwc_simulations_results -  select_d= ", " select_n= ", " select_rho= ", " select_sigma= ", " select_chain= ", " seed_sigma= ", " seed_data= ", ".Rdata")

pattern <- paste(str_escape(word_to_remove), collapse = "|")
ff_cw_clean <- str_remove_all(ff_cw, pattern)
split_list <- strsplit(ff_cw_clean, " ")
mat_cw <- do.call(rbind, split_list)
data_cw <- as.data.frame(mat_cw)
colnames(data_cw) <- c("do_best_init", "select_d", "select_n", "select_kappa", "select_sigma", "select_chain", "seed_sigma", "seed_data")
i <- 1
data_cw[[i]] <- as.logical(as.character(data_cw[[i]]))
for (i in 2:ncol(data_cw)) {
  data_cw[[i]] <- as.numeric(as.character(data_cw[[i]]))
}
data_cw[, 1] <- FALSE

# Containers for diagnostics and true parameters.
list_ret_cw <- list()
list_ret_pn <- list()


h <- 0


old_wd <- getwd()
setwd(paste(old_wd, "simulations/output", sep = "/"))
dir_out <- paste(old_wd, "simulations/output_diagnostic/", sep = "/")
paramaters_cw <- list()
paramaters_pn <- list()


# Loop over simulation scenarios. Within each scenario the script loads all
# available chains, computes diagnostics for the copula model, then repeats the
# same calculations for the TPN model.
for (init_sel in c(FALSE))
{
  for (d_sel in c(4,3, 2, 1))
  {
    for (n_sel in c(1, 2))
    {
      for (kappa_sel in 4:4)
      {
        rho_sel <- kappa_sel
        for (select_sigma in 1:2)
        {
          for (seed_data in 1:25)
          {
            h <- h + 1
            print(paste("do_best_init=", init_sel, "select_d=", d_sel, "select_n=", n_sel, "select_kappa=", kappa_sel, "select_sigma=", select_sigma, "seed_data=", seed_data))


            w_data <- which(data_cw$do_best_init == init_sel & data_cw$select_d == d_sel & data_cw$select_n == n_sel & data_cw$select_kappa == kappa_sel & select_sigma == data_cw$select_sigma & data_cw$seed_data == seed_data)
            data_plot <- data_cw[w_data, ]

            ff_plot <- ff_cw[w_data]


            list_ret_cw[[h]] <- list()

            list_ret_cw[[h]]$parameters <- list(init_sel = init_sel, d_sel = d_sel, n_sel = n_sel, kappa_sel = kappa_sel, select_sigma = select_sigma, seed_data = seed_data)
            if (length(ff_plot) < 2) {
              print("QUI CW")
              print(ff_plot)
              print(paste("do_best_init=", init_sel, "select_d=", d_sel, "select_n=", n_sel, "select_rho=", kappa_sel, "select_sigma=", select_sigma, "seed_data=", seed_data))
            } else {
              # Load all available copula-model chains for this scenario.
              chain <- list()
              load(ff_plot[1])
              chain[[1]] <- res_list
              load(ff_plot[2])
              chain[[2]] <- res_list

              if (length(ff_plot) >= 3) {
                load(ff_plot[3])
                chain[[3]] <- res_list
              }
              if (length(ff_plot) >= 4) {
                load(ff_plot[4])
                chain[[4]] <- res_list
              }
              if (length(ff_plot) >= 5) {
                load(ff_plot[5])
                chain[[5]] <- res_list
              }

              d <- chain[[1]]$d

              # Compute diagnostics for the unique off-diagonal covariance
              # elements.
              nsim <- nrow(chain[[1]]$sigma_s_out)

              index_non_1 <- which(lower.tri(matrix(0, d, d)))
              paramaters_cw[[h]] <- list(Sigma_s = res_list$Sigma_s, Sigma_c = res_list$Sigma_c, rho = res_list$rho, mu = res_list$mu, pos_def = res_list$pos_def, index_non_1 = index_non_1)
              chains <- lapply(chain, function(ch) ch[["sigma_s_out"]])

              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )


              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              rhsat_sigma <- cbind(rhat_vec, 1)
              list_ret_cw[[h]]$rhat_sigma <- cbind(rhat_vec[index_non_1], 1)
              list_ret_cw[[h]]$ess_sigma <- ess_vec[index_non_1]
              list_ret_cw[[h]]$ess_bulk_sigma <- ess_bulk_vec[index_non_1]
              list_ret_cw[[h]]$ess_tail_sigma <- ess_tail_vec[index_non_1]


              # posterior::rhat(draws)

              # ess_bulk(draws)

              # ess_tail(draws)
              # rhat <- gelman.diag(mcmc_list, multivariate = FALSE)

              # rhsat_sigma <- rhat$psrf


              # ess <- effectiveSize(mcmc_list)

              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$sigma_s_out)
              )
              list_ret_cw[[h]]$multiESS_sigma <- multiESS_base2(sim_par[, index_non_1], max_lag = 100)$mESS
              list_ret_cw[[h]]$multiESS_sigma_stand <- list_ret_cw[[h]]$multiESS_sigma / ((d^2 - d) / 2)
              # sim_par <- rbind(chain1$sigma_s_out, chain2$sigma_s_out)

              true_par <- matrix(c(chain[[1]]$Sigma_s), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- sim_par - true_par
              # R x d matrix
              bias <- colMeans(sim_par - true_par)
              mse <- colMeans((sim_par - true_par)^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_cw[[h]]$bias_sigma <- bias[index_non_1]
              list_ret_cw[[h]]$mse_sigma <- mse[index_non_1]
              list_ret_cw[[h]]$coverage_sigma <- coverage[index_non_1]
              list_ret_cw[[h]]$ICL_sigma <- qq3 - qq1

              ## SIGMA C
              # ! SIGMA
              # x <- abind::abind(
              #  chain1$sigma_c_out,
              #  chain2$sigma_c_out,
              #  along = 3
              # )
              # x <- aperm(x, c(1, 3, 2))
              # dimnames(x) <- list(
              #  NULL,
              #  paste0("chain", 1:2),
              #  paste0("par", 1:dim(x)[3])
              # )

              chains <- lapply(chain, function(ch) ch[["sigma_c_out"]])

              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )
              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              rhsat_sigma_c <- cbind(rhat_vec, 1)
              list_ret_cw[[h]]$rhat_sigma_c <- cbind(rhat_vec[index_non_1], 1)
              list_ret_cw[[h]]$ess_sigma_c <- ess_vec[index_non_1]
              list_ret_cw[[h]]$ess_bulk_sigma_c <- ess_bulk_vec[index_non_1]
              list_ret_cw[[h]]$ess_tail_sigma_c <- ess_tail_vec[index_non_1]
              # mcmc_list <- mcmc.list(
              #  mcmc(chain1$sigma_c_out),
              #  mcmc(chain2$sigma_c_out)
              # )
              # rhat <- gelman.diag(mcmc_list, multivariate = FALSE)
              # list_ret_cw[[h]]$rhat_sigma_c <- rhat$psrf[index_non_1,]
              # rhsat_sigma_c <- rhat$psrf


              # ess <- effectiveSize(mcmc_list)

              # list_ret_cw[[h]]$ess_sigma_c <- ess[index_non_1]

              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$sigma_c_out)
              )
              list_ret_cw[[h]]$multiESS_sigma_c <- multiESS_base2(sim_par[, index_non_1], max_lag = 100)$mESS
              list_ret_cw[[h]]$multiESS_sigma_c_stand <- list_ret_cw[[h]]$multiESS_sigma_c / ((d^2 - d) / 2)

              true_par <- matrix(c(chain[[1]]$Sigma_c), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- sim_par - true_par
              # R x d matrix
              bias <- colMeans(sim_par - true_par)
              mse <- colMeans((sim_par - true_par)^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_cw[[h]]$bias_sigma_c <- bias[index_non_1]
              list_ret_cw[[h]]$mse_sigma_c <- mse[index_non_1]
              list_ret_cw[[h]]$coverage_sigma_c <- coverage[index_non_1]
              list_ret_cw[[h]]$ICL_sigma_c <- qq3 - qq1


              # ! rho
              # x <- abind::abind(
              #  chain1$rho_out,
              #  chain2$rho_out,
              #  along = 3
              # )
              # x <- aperm(x, c(1, 3, 2))
              # dimnames(x) <- list(
              #  NULL,
              #  paste0("chain", 1:2),
              #  paste0("par", 1:dim(x)[3])
              # )
              chains <- lapply(chain, function(ch) ch[["rho_out"]])

              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )

              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              list_ret_cw[[h]]$rhat_rho <- cbind(rhat_vec, 1)
              list_ret_cw[[h]]$ess_rho <- ess_vec
              list_ret_cw[[h]]$ess_bulk_rho <- ess_bulk_vec
              list_ret_cw[[h]]$ess_tail_rho <- ess_tail_vec

              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$rho_out)
              )

              true_par <- matrix(c(chain[[1]]$rho), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- sim_par - true_par
              # R x d matrix
              bias <- colMeans(sim_par - true_par)
              mse <- colMeans((sim_par - true_par)^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_cw[[h]]$bias_rho <- bias
              list_ret_cw[[h]]$mse_rho <- mse
              list_ret_cw[[h]]$coverage_rho <- coverage
              list_ret_cw[[h]]$ICL_rho <- qq3 - qq1

              # ! mu
              mu <- chain[[1]]$mu # nsim <- nrow(chain1$mu_out)
              # tot_mean <- rbind(chain1$mu_out, chain2$mu_out)
              tot_mean <- do.call(
                rbind,
                lapply(chain, function(ch) ch$mu_out)
              )
              circ_mean_vec <- c()
              for (j in 1:res_list$d)
              {
                circ_mean <- atan2(mean(sin(tot_mean[, j])), mean(cos(tot_mean[, j])))
                tot_mean[, j] <- (tot_mean[, j] - circ_mean + pi) %% (2 * pi)
                circ_mean_vec[j] <- circ_mean
              }
              for (iii in 1:length(chain))
              {
                chain[[iii]]$mu_app <- tot_mean[((iii - 1) * nsim + 1):(iii * nsim), ]
              }
              # chain1$mu_app <- tot_mean[1:nsim, ]
              # chain2$mu_app <- tot_mean[(nsim + 1):(2 * nsim), ]


              # mcmc_list <- mcmc.list(
              #  mcmc(chain1$mu_app),
              #  mcmc(chain2$mu_app)
              # )
              # rhat <- gelman.diag(mcmc_list, multivariate = FALSE)
              # list_ret_cw[[h]]$rhat_mu <- rhat$psrf


              # ess <- effectiveSize(mcmc_list)

              # list_ret_cw[[h]]$ess_mu <- ess
              chains <- lapply(chain, function(ch) ch[["mu_app"]])

              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )

              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              list_ret_cw[[h]]$rhat_mu <- cbind(rhat_vec, 1)
              list_ret_cw[[h]]$ess_mu <- ess_vec
              list_ret_cw[[h]]$ess_bulk_mu <- ess_bulk_vec
              list_ret_cw[[h]]$ess_tail_mu <- ess_tail_vec

              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$mu_app)
              )

              true_par <- matrix(c((mu - circ_mean_vec + pi) %% (2 * pi)), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- atan2(sin(sim_par - true_par), cos(sim_par - true_par))
              # R x d matrix
              bias <- colMeans(centered_par)
              mse <- colMeans(centered_par^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_cw[[h]]$bias_mu <- bias
              list_ret_cw[[h]]$mse_mu <- mse
              list_ret_cw[[h]]$coverage_mu <- coverage
              list_ret_cw[[h]]$ICL_mu <- qq3 - qq1

              ## ! PLOTS
              # pdf(paste(dir_out, name, "new_cw_post_", paste("do_best_init=", init_sel, "select_d=", d_sel, "select_n=", n_sel, "select_kappa=", kappa_sel, "select_sigma=", select_sigma, "seed_data=", seed_data), ".pdf", sep = ""))
              # par(mfrow = c(3, 3))
              # for (id in 1:d)
              # {
              #  plot(chain[[1]]$mu_app[, id], type = "l", main = round(chain[[1]]$mu[id], 3))
              #  for (iii in 2:length(chain))
              #  {
              #    lines(chain[[iii]]$mu_app[, id], col = iii)
              #  }
              #  abline(h = pi, col = 2)
              #  plot(density(chain[[1]]$mu_app[, id]))
              #  for (iii in 2:length(chain))
              #  {
              #    lines(density(chain[[iii]]$mu_app[, id]), col = iii)
              #  }
              #  abline(v = pi, col = 2)
              #  plot_multi_acf(chain, "mu_app", id, max_lag = 20)
              # }
              # par(mfrow = c(3, 3))
              # for (id in 1:d) {
              #  plot(chain[[1]]$rho_out[, id], type = "l", main = round(chain[[1]]$rho[id], 3))
              #  for (iii in 2:length(chain)) {
              #    lines(chain[[iii]]$rho_out[, id], col = iii)
              #  }
              #  abline(h = chain[[1]]$rho[id], col = 2)

              #  plot(density(chain[[1]]$rho_out[, id]))
              #  for (iii in 2:length(chain)) {
              #    lines(density(chain[[iii]]$rho_out[, id]), col = iii)
              #  }
              #  abline(v = chain[[1]]$rho[id], col = 2)

              #  plot_multi_acf(chain, "rho_out", id, max_lag = 20)
              # }

              # par(mfrow = c(3, 3))
              # hhhh <- 1
              # for (id in 1:d) {
              #  for (jd in 1:d) {
              #    if (jd < id) {
              #      plot(chain[[1]]$sigma_s_out[, hhhh],
              #        type = "l",
              #        main = round(chain[[1]]$Sigma_s[id, jd], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(chain[[iii]]$sigma_s_out[, hhhh], col = iii)
              #      }
              #      abline(h = chain[[1]]$Sigma_s[id, jd], col = 2)

              #      plot(density(chain[[1]]$sigma_s_out[, hhhh]),
              #        main = round(rhsat_sigma[hhhh, 1], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(density(chain[[iii]]$sigma_s_out[, hhhh]), col = iii)
              #      }
              #      abline(v = chain[[1]]$Sigma_s[id, jd], col = 2)

              #      plot_multi_acf(chain, "sigma_s_out", hhhh, max_lag = 20, par_title = paste("Sigma_s[", id, ",", jd, "]"))
              #    }


              #    hhhh <- hhhh + 1
              #  }
              # }

              # par(mfrow = c(3, 3))
              # hhhh <- 1
              # for (id in 1:d) {
              #  for (jd in 1:d) {
              #    if (jd < id) {
              #      plot(chain[[1]]$sigma_c_out[, hhhh],
              #        type = "l",
              #        main = round(chain[[1]]$Sigma_c[id, jd], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(chain[[iii]]$sigma_c_out[, hhhh], col = iii)
              #      }
              #      abline(h = chain[[1]]$Sigma_c[id, jd], col = 2)

              #      plot(density(chain[[1]]$sigma_c_out[, hhhh]),
              #        main = round(rhsat_sigma_c[hhhh, 1], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(density(chain[[iii]]$sigma_c_out[, hhhh]), col = iii)
              #      }
              #      abline(v = chain[[1]]$Sigma_c[id, jd], col = 2)

              #      plot_multi_acf(chain, "sigma_c_out", hhhh, max_lag = 20, par_title = paste("Sigma_c[", id, ",", jd, "]"))
              #    }


              #    hhhh <- hhhh + 1
              #  }
              # }
              # dev.off()
            }

            # Repeat the same diagnostics for the TPN model under the matching
            # simulation scenario.
            w_data <- which(data_pn$do_best_init == init_sel & data_pn$select_d == d_sel & data_pn$select_n == n_sel & data_pn$select_kappa == kappa_sel & select_sigma == data_pn$select_sigma & data_pn$seed_data == seed_data)
            data_plot <- data_pn[w_data, ]

            ff_plot <- ff_pn[w_data]


            list_ret_pn[[h]] <- list()
            list_ret_pn[[h]]$parameters <- list(init_sel = init_sel, d_sel = d_sel, n_sel = n_sel, kappa_sel = kappa_sel, select_sigma = select_sigma, seed_data = seed_data)
            if (length(ff_plot) < 2) {
              print("QUI PN")
              print(ff_plot)
              print(paste("do_best_init=", init_sel, "select_d=", d_sel, "select_n=", n_sel, "select_kappa=", kappa_sel, "select_sigma=", select_sigma, "seed_data=", seed_data))
            } else {
              # Load all available TPN chains for this scenario.
              chain <- list()
              load(ff_plot[1])
              chain[[1]] <- res_list
              load(ff_plot[2])
              chain[[2]] <- res_list

              if (length(ff_plot) >= 3) {
                load(ff_plot[3])
                chain[[3]] <- res_list
              }
              if (length(ff_plot) >= 4) {
                load(ff_plot[4])
                chain[[4]] <- res_list
              }
              if (length(ff_plot) >= 5) {
                load(ff_plot[5])
                chain[[5]] <- res_list
              }

              d <- chain[[1]]$d

              ## chain 1
              nsim <- nrow(chain[[1]]$sigma_s_out)

              index_non_1 <- which(lower.tri(matrix(0, d, d)))
              paramaters_pn[[h]] <- list(Sigma_s = res_list$Sigma_s, Sigma_c = res_list$Sigma_c, kappa = res_list$kappa, mu = res_list$mu, pos_def = res_list$pos_def, index_non_1 = index_non_1)
              chains <- lapply(chain, function(ch) ch[["sigma_s_out"]])

              # ! SIGMA
              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )

              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              rhsat_sigma <- cbind(rhat_vec, 1)
              list_ret_pn[[h]]$rhat_sigma <- cbind(rhat_vec[index_non_1], 1)
              list_ret_pn[[h]]$ess_sigma <- ess_vec[index_non_1]
              list_ret_pn[[h]]$ess_bulk_sigma <- ess_bulk_vec[index_non_1]
              list_ret_pn[[h]]$ess_tail_sigma <- ess_tail_vec[index_non_1]


              # posterior::rhat(draws)

              # ess_bulk(draws)

              # ess_tail(draws)
              # rhat <- gelman.diag(mcmc_list, multivariate = FALSE)

              # rhsat_sigma <- rhat$psrf


              # ess <- effectiveSize(mcmc_list)


              # sim_par <- rbind(chain1$sigma_s_out, chain2$sigma_s_out)
              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$sigma_s_out)
              )
              list_ret_pn[[h]]$multiESS_sigma <- multiESS_base2(sim_par[, index_non_1], max_lag = 100)$mESS
              list_ret_pn[[h]]$multiESS_sigma_stand <- list_ret_pn[[h]]$multiESS_sigma / ((d^2 - d) / 2)

              true_par <- matrix(c(chain[[1]]$Sigma_s), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- sim_par - true_par
              # R x d matrix
              bias <- colMeans(sim_par - true_par)
              mse <- colMeans((sim_par - true_par)^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_pn[[h]]$bias_sigma <- bias[index_non_1]
              list_ret_pn[[h]]$mse_sigma <- mse[index_non_1]
              list_ret_pn[[h]]$coverage_sigma <- coverage[index_non_1]
              list_ret_pn[[h]]$ICL_sigma <- qq3 - qq1

              ## SIGMA C
              # ! SIGMA
              chains <- lapply(chain, function(ch) ch[["sigma_c_out"]])

              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )

              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              rhsat_sigma_c <- cbind(rhat_vec, 1)
              list_ret_pn[[h]]$rhat_sigma_c <- cbind(rhat_vec[index_non_1], 1)
              list_ret_pn[[h]]$ess_sigma_c <- ess_vec[index_non_1]
              list_ret_pn[[h]]$ess_bulk_sigma_c <- ess_bulk_vec[index_non_1]
              list_ret_pn[[h]]$ess_tail_sigma_c <- ess_tail_vec[index_non_1]
              # mcmc_list <- mcmc.list(
              #  mcmc(chain1$sigma_c_out),
              #  mcmc(chain2$sigma_c_out)
              # )
              # rhat <- gelman.diag(mcmc_list, multivariate = FALSE)
              # list_ret_cw[[h]]$rhat_sigma_c <- rhat$psrf[index_non_1,]
              # rhsat_sigma_c <- rhat$psrf


              # ess <- effectiveSize(mcmc_list)

              # list_ret_cw[[h]]$ess_sigma_c <- ess[index_non_1]

              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$sigma_c_out)
              )
              list_ret_pn[[h]]$multiESS_sigma_c <- multiESS_base2(sim_par[, index_non_1], max_lag = 100)$mESS
              list_ret_pn[[h]]$multiESS_sigma_c_stand <- list_ret_pn[[h]]$multiESS_sigma_c / ((d^2 - d) / 2)


              true_par <- matrix(c(chain[[1]]$Sigma_c), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- sim_par - true_par
              # R x d matrix
              bias <- colMeans(sim_par - true_par)
              mse <- colMeans((sim_par - true_par)^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_pn[[h]]$bias_sigma_c <- bias[index_non_1]
              list_ret_pn[[h]]$mse_sigma_c <- mse[index_non_1]
              list_ret_pn[[h]]$coverage_sigma_c <- coverage[index_non_1]
              list_ret_pn[[h]]$ICL_sigma_c <- qq3 - qq1


              # ! rho
              chains <- lapply(chain, function(ch) ch[["kappa_out"]])

              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )

              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              list_ret_pn[[h]]$rhat_kappa <- cbind(rhat_vec, 1)
              list_ret_pn[[h]]$ess_kappa <- ess_vec
              list_ret_pn[[h]]$ess_bulk_kappa <- ess_bulk_vec
              list_ret_pn[[h]]$ess_tail_kappa <- ess_tail_vec

              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$kappa_out)
              )

              true_par <- matrix(c(chain[[1]]$kappa), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- sim_par - true_par
              # R x d matrix
              bias <- colMeans(sim_par - true_par)
              mse <- colMeans((sim_par - true_par)^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_pn[[h]]$bias_kappa <- bias
              list_ret_pn[[h]]$mse_kappa <- mse
              list_ret_pn[[h]]$coverage_kappa <- coverage
              list_ret_pn[[h]]$ICL_kappa <- qq3 - qq1

              # ! mu
              mu <- chain[[1]]$mu # nsim <- nrow(chain1$mu_out)
              # tot_mean <- rbind(chain1$mu_out, chain2$mu_out)
              tot_mean <- do.call(
                rbind,
                lapply(chain, function(ch) ch$mu_out)
              )
              circ_mean_vec <- c()
              for (j in 1:res_list$d)
              {
                circ_mean <- atan2(mean(sin(tot_mean[, j])), mean(cos(tot_mean[, j])))
                tot_mean[, j] <- (tot_mean[, j] - circ_mean + pi) %% (2 * pi)
                circ_mean_vec[j] <- circ_mean
              }
              for (iii in 1:length(chain))
              {
                chain[[iii]]$mu_app <- tot_mean[((iii - 1) * nsim + 1):(iii * nsim), ]
              }


              # mcmc_list <- mcmc.list(
              #  mcmc(chain1$mu_app),
              #  mcmc(chain2$mu_app)
              # )
              # rhat <- gelman.diag(mcmc_list, multivariate = FALSE)
              # list_ret_cw[[h]]$rhat_mu <- rhat$psrf


              # ess <- effectiveSize(mcmc_list)

              # list_ret_cw[[h]]$ess_mu <- ess
              chains <- lapply(chain, function(ch) ch[["mu_app"]])

              x <- do.call(abind::abind, c(chains, along = 3))
              x <- aperm(x, c(1, 3, 2))
              dimnames(x) <- list(
                NULL,
                paste0("chain", 1:length(chain)),
                paste0("par", 1:dim(x)[3])
              )

              draws <- as_draws_array(x)

              rhat_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)

                posterior::rhat(mat)
              })
              ess_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_basic(mat)
              })
              ess_bulk_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_bulk(mat)
              })
              ess_tail_vec <- sapply(posterior::variables(draws), function(v) {
                mat <- posterior::extract_variable_matrix(draws, variable = v)
                posterior::ess_tail(mat)
              })

              list_ret_pn[[h]]$rhat_mu <- cbind(rhat_vec, 1)
              list_ret_pn[[h]]$ess_mu <- ess_vec
              list_ret_pn[[h]]$ess_bulk_mu <- ess_bulk_vec
              list_ret_pn[[h]]$ess_tail_mu <- ess_tail_vec

              sim_par <- do.call(
                rbind,
                lapply(chain, function(ch) ch$mu_app)
              )

              true_par <- matrix(c((mu - circ_mean_vec + pi) %% (2 * pi)), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
              centered_par <- atan2(sin(sim_par - true_par), cos(sim_par - true_par))
              # R x d matrix
              bias <- colMeans(centered_par)
              mse <- colMeans(centered_par^2)

              qq1 <- apply(centered_par, 2, quantile, probs = 0.025)
              qq3 <- apply(centered_par, 2, quantile, probs = 0.975)
              coverage <- (qq1 <= 0) & (qq3 >= 0)

              list_ret_pn[[h]]$bias_mu <- bias
              list_ret_pn[[h]]$mse_mu <- mse
              list_ret_pn[[h]]$coverage_mu <- coverage
              list_ret_pn[[h]]$ICL_mu <- qq3 - qq1

              ## ! PLOTS
              # pdf(paste(dir_out, name, "new_pn_post_", paste("do_best_init=", init_sel, "select_d=", d_sel, "select_n=", n_sel, "select_kappa=", kappa_sel, "select_sigma=", select_sigma, "seed_data=", seed_data), ".pdf", sep = ""))
              # par(mfrow = c(3, 3))
              # for (id in 1:d)
              # {
              #  plot(chain[[1]]$mu_app[, id], type = "l", main = round(chain[[1]]$mu[id], 3))
              #  for (iii in 2:length(chain))
              #  {
              #    lines(chain[[iii]]$mu_app[, id], col = iii)
              #  }
              #  abline(h = pi, col = 2)
              #  plot(density(chain[[1]]$mu_app[, id]))
              #  for (iii in 2:length(chain))
              #  {
              #    lines(density(chain[[iii]]$mu_app[, id]), col = iii)
              #  }
              #  abline(v = pi, col = 2)
              #  plot_multi_acf(chain, "mu_app", id, max_lag = 20)
              # }
              # par(mfrow = c(3, 3))
              # for (id in 1:d) {
              #  plot(chain[[1]]$kappa_out[, id], type = "l", main = round(chain[[1]]$kappa[id], 3))
              #  for (iii in 2:length(chain)) {
              #    lines(chain[[iii]]$kappa_out[, id], col = iii)
              #  }
              #  abline(h = chain[[1]]$kappa[id], col = 2)

              #  plot(density(chain[[1]]$kappa_out[, id]))
              #  for (iii in 2:length(chain)) {
              #    lines(density(chain[[iii]]$kappa_out[, id]), col = iii)
              #  }
              #  abline(v = chain[[1]]$kappa[id], col = 2)

              #  plot_multi_acf(chain, "kappa_out", id, max_lag = 20)
              # }

              # par(mfrow = c(3, 3))
              # hhhh <- 1
              # for (id in 1:d) {
              #  for (jd in 1:d) {
              #    if (jd < id) {
              #      plot(chain[[1]]$sigma_s_out[, hhhh],
              #        type = "l",
              #        main = round(chain[[1]]$Sigma_s[id, jd], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(chain[[iii]]$sigma_s_out[, hhhh], col = iii)
              #      }
              #      abline(h = chain[[1]]$Sigma_s[id, jd], col = 2)

              #      plot(density(chain[[1]]$sigma_s_out[, hhhh]),
              #        main = round(rhsat_sigma[hhhh, 1], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(density(chain[[iii]]$sigma_s_out[, hhhh]), col = iii)
              #      }
              #      abline(v = chain[[1]]$Sigma_s[id, jd], col = 2)

              #      plot_multi_acf(chain, "sigma_s_out", hhhh, max_lag = 20)
              #    }


              #    hhhh <- hhhh + 1
              #  }
              # }

              # par(mfrow = c(3, 3))
              # hhhh <- 1
              # for (id in 1:d) {
              #  for (jd in 1:d) {
              #    if (jd < id) {
              #      plot(chain[[1]]$sigma_c_out[, hhhh],
              #        type = "l",
              #        main = round(chain[[1]]$Sigma_c[id, jd], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(chain[[iii]]$sigma_c_out[, hhhh], col = iii)
              #      }
              #      abline(h = chain[[1]]$Sigma_c[id, jd], col = 2)

              #      plot(density(chain[[1]]$sigma_c_out[, hhhh]),
              #        main = round(rhsat_sigma_c[hhhh, 1], 3)
              #      )
              #      for (iii in 2:length(chain)) {
              #        lines(density(chain[[iii]]$sigma_c_out[, hhhh]), col = iii)
              #      }
              #      abline(v = chain[[1]]$Sigma_c[id, jd], col = 2)

              #      plot_multi_acf(chain, "sigma_c_out", hhhh, max_lag = 20)
              #    }


              #    hhhh <- hhhh + 1
              #  }
              #}
              #dev.off()
            }
            # Save the diagnostic summaries for all scenarios processed so far.
            save(list_ret_cw, list_ret_pn, paramaters_cw, paramaters_pn, file = paste(dir_out, name, "new_post_analisi_results_stat.Rdata", sep = ""))
          }
        }
      }
    }
  }
}

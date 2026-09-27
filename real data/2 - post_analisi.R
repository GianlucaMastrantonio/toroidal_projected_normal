tex_quantile <- function(quantile, x) {
  paste0(
    "(",
    paste(
      round(quantile(x,
        prob = quantile
      ), 3),
      collapse = " "
    ),
    ")"
  )
}


library(glue)
# tables
# Fill these matrices
Sigma <- matrix("", 3, 4)
Sigma_c <- matrix("", 3, 4)
Mu <- matrix("", 3, 4)
Kappa <- matrix("", 3, 4)

# Example:
# Sigma[1,] <- c("1.01","1.00","1.00","1.00","1.00","1.00")

row_tex <- function(label, metric, vals) {
  paste0(
    label, "&", metric, "&",
    paste(vals, collapse = " & "),
    "\\\\"
  )
}

# row_tex <- function(label, metric, vals) {
#  vals <- unname(as.character(vals))

#  paste0(label, " & ", metric, " & ", paste(vals, collapse = " & "), "\\\\")
# }


table_tex <- function(Sigma,
                      Sigma_c,
                      Mu,
                      Kappa,
                      groups = c("Independent Model", "Full Model"),
                      n_vals = c(35, 75, 35, 75),
                      caption = "PN sigma dipendente.", prec_name = "\\boldsymbol{\\kappa}") {
  tex <- c(
    "\\begin{table}[t]",
    "\\centering",
    "\\scriptsize",
    "\\begin{tabular}{c|c|cc|cc}",
    "\\hline",
    paste0(
      "&   & ",
      "\\multicolumn{2}{c|}{", groups[1], "} & ",
      "\\multicolumn{2}{c}{", groups[2], "}\\\\"
    ),
    "\\cline{1-6}",
    paste0(
      "& d & ",
      paste(sprintf("$%s$", n_vals), collapse = " & "),
      " \\\\"
    ),
    "\\hline",
    row_tex("", "$\\hat{R}$", Sigma[1, ]),
    row_tex("$\\boldsymbol{\\Sigma}$", "ESS", Sigma[2, ]),
    # row_tex("$\\boldsymbol{\\Sigma}$", "Coverage", Sigma[3, ]),
    # row_tex("", "Bias", Sigma[4, ]),
    # row_tex("", "MSE", Sigma[5, ]),
    row_tex("", "CIL", Sigma[3, ]),
    "\\hline",
    row_tex("", "$\\hat{R}$", Sigma_c[1, ]),
    row_tex("$\\boldsymbol{\\Sigma}_c$", "ESS", Sigma_c[2, ]),
    # row_tex("$\\boldsymbol{\\Sigma}_c$", "Coverage", Sigma_c[3, ]),
    # row_tex("", "Bias", Sigma_c[4, ]),
    # row_tex("", "MSE", Sigma_c[5, ]),
    row_tex("", "CIL", Sigma_c[3, ]),
    "\\hline",
    row_tex("", "$\\hat{R}$", Mu[1, ]),
    row_tex("$\\boldsymbol{\\mu}$", "ESS", Mu[2, ]),
    # row_tex("$\\boldsymbol{\\mu}$", "Coverage", Mu[3, ]),
    # row_tex("", "Bias", Mu[4, ]),
    # row_tex("", "MSE", Mu[5, ]),
    row_tex("", "CIL", Mu[3, ]),
    "\\hline",
    row_tex("", "$\\hat{R}$", Kappa[1, ]),
    row_tex(paste0("$", prec_name, "$"), "ESS", Kappa[2, ]),
    # row_tex(paste0("$", prec_name, "$"), "Coverage", Kappa[3, ]),
    # row_tex("", "Bias", Kappa[4, ]),
    # row_tex("", "MSE", Kappa[5, ]),
    row_tex("", "CIL", Kappa[3, ]),
    "\\hline",
    "\\end{tabular}",
    paste0("\\caption{", caption, "}"),
    "\\end{table}"
  )

  paste(tex, collapse = "\n")
}


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
    main = if (is.null(par_title)) paste(par_name, "par", par_id) else par_title
  )

  for (cc in 2:length(acf_list)) {
    lines(acf_list[[cc]]$lag, acf_list[[cc]]$acf, col = cc)
  }

  abline(h = 0, lty = 2)
}
library(stringr)
library(coda)
library(posterior)

library(stringr)

dir_data <- "real data/output/"
dir_out <- "real data/output_diagnostic/"
ff <- list.files(dir_data)



name_sim <- "REAL"
w <- grep(".Rdata", ff)
ff <- ff[w]
w <- grep(name_sim, ff)
ff <- ff[w]


df <- data.frame(
  #  ind_ = str_extract(ff, "(?<=IND)[A-Z]+_[A-Z]+(?=[0-9])"),
  ind = str_extract(ff, "(?<=IND)(TRUE|FALSE)(?=_)"),
  only_ess = str_extract(ff, "(?<=IND(?:TRUE|FALSE)_)(TRUE|FALSE)(?=[0-9])"),
  ntry = as.numeric(str_extract(ff, "(?<=IND(?:TRUE|FALSE)_(?:TRUE|FALSE))[0-9]+")) %% (100),
  type_ess = (as.numeric(str_extract(ff, "(?<=IND(?:TRUE|FALSE)_(?:TRUE|FALSE))[0-9]+")) - as.numeric(str_extract(ff, "(?<=IND(?:TRUE|FALSE)_(?:TRUE|FALSE))[0-9]+")) %% (100)) / 100,
  do_best_init = str_extract(ff, "(?<=do_best_init=)(TRUE|FALSE)") == "TRUE",
  do_small = as.numeric(str_extract(ff, "(?<=do_small=T)[0-9]")),
  molt_iter = str_extract(ff, "(?<=molt_iter=)(1|2)"),
  model = str_extract(ff, "(cwc|tpn)(?=_seed)"),
  seed = as.numeric(str_extract(ff, "(?<=_seed)[0-9]+"))
)


df$setting_id <- as.integer(
  interaction(df[, c("ind", "only_ess", "ntry", "do_small", "molt_iter", "model", "do_best_init")],
    drop = TRUE
  )
)

crps_circ_old <- function(y, x) {
  d1 <- 1 - cos(y - x)

  dx <- outer(x, x, function(a, b) 1 - cos(a - b))

  mean(d1) - 0.5 * mean(dx)
}
crps_circ <- function(y, x) {
  cx <- cos(x)
  sx <- sin(x)
  L <- length(x)

  mean_d1 <- 1 - cos(y) * mean(cx) - sin(y) * mean(sx)

  R2 <- sum(cx)^2 + sum(sx)^2
  mean_dx <- 1 - R2 / L^2

  mean_d1 - 0.5 * mean_dx
}

# df[df$only_ess == TRUE & df$do_small == 25 & df$model == "tpn", ]

un_id <- unique(df$setting_id)
un_id <- un_id[order(un_id)]
###
list_ret <- list()
set_id_unique <- c()

set.seed(1)
# load("/Users/gianlucamastrantonio/Politecnico di Torino Staff Dropbox/Gianluca Mastrantonio/lavori/gitrepo/toroidal_projected_normal/real data/data/data_stations_code.RData")

# new_n_miss <- floor(nrow(theta) * 0.1 * ncol(theta))
# new_set_miss <- data.frame(obs = 1:nrow(theta), d = rep(1:ncol(theta), each = nrow(theta)))
# new_index_miss <- sample(1:nrow(new_set_miss), new_n_miss, replace = FALSE)
# new_index_miss <- new_index_miss[order(new_index_miss)]
# new_set_miss <- new_set_miss[new_index_miss, ]


# new_theta_miss <- c()
# for (ii in 1:nrow(new_set_miss))
# {
#  iobs <- new_set_miss[ii, 1]
#  id <- new_set_miss[ii, 2]
#  new_theta_miss[ii] <- theta[iobs, id]
# }



for (id_real in un_id)
{
  w_which_df <- df$setting_id == id_real
  set_id_code <- 1:5


  ex_df <- df[which(w_which_df)[1], ]
  ff_plot <- ff[which(df$setting_id == id_real)][set_id_code]
  set_id_unique[id_real] <- c(which(w_which_df)[1])
  print(id_real)
  print(ex_df)

  chain <- list()


  if (length(ff_plot) > 1) {
    load(paste(dir_data, ff_plot[1], sep = ""))
    if (ex_df$model[1] == "tpn") {
      prec <- kappa_out
    } else {
      prec <- rho_out
    }
    res_list <- list(
      circ_mean = mu_out,
      prec = prec,
      sigma_s_out = sigma_s_out,
      sigma_c_out = sigma_c_out,
      missig_out = out_mcmc$missig_out
    )
    chain[[1]] <- res_list
    load(paste(dir_data, ff_plot[2], sep = ""))

    if (ex_df$model[1] == "tpn") {
      prec <- kappa_out
    } else {
      prec <- rho_out
    }
    res_list <- list(
      circ_mean = mu_out,
      prec = prec,
      sigma_s_out = sigma_s_out,
      sigma_c_out = sigma_c_out,
      missig_out = out_mcmc$missig_out
    )
    chain[[2]] <- res_list

    if (length(ff_plot) >= 3) {
      load(paste(dir_data, ff_plot[3], sep = ""))

      if (ex_df$model[1] == "tpn") {
        prec <- kappa_out
      } else {
        prec <- rho_out
      }
      res_list <- list(
        circ_mean = mu_out,
        prec = prec,
        sigma_s_out = sigma_s_out,
        sigma_c_out = sigma_c_out,
        missig_out = out_mcmc$missig_out
      )
      chain[[3]] <- res_list
    }
    if (length(ff_plot) >= 4) {
      load(paste(dir_data, ff_plot[4], sep = ""))

      if (ex_df$model[1] == "tpn") {
        prec <- kappa_out
      } else {
        prec <- rho_out
      }
      res_list <- list(
        circ_mean = mu_out,
        prec = prec,
        sigma_s_out = sigma_s_out,
        sigma_c_out = sigma_c_out,
        missig_out = out_mcmc$missig_out
      )
      chain[[4]] <- res_list
    }
    if (length(ff_plot) >= 5) {
      load(paste(dir_data, ff_plot[5], sep = ""))

      if (ex_df$model[1] == "tpn") {
        prec <- kappa_out
      } else {
        prec <- rho_out
      }
      res_list <- list(
        circ_mean = mu_out,
        prec = prec,
        sigma_s_out = sigma_s_out,
        sigma_c_out = sigma_c_out,
        missig_out = out_mcmc$missig_out
      )
      chain[[5]] <- res_list
    }
    ### CRPS
    # index_miss <- sample(1:n, n_miss, replace = FALSE)
    # index_miss <- index_miss[order(index_miss)]
    # theta_miss <- theta_all[index_miss, ]
    # theta <- theta_all[-index_miss, ]
    #source("/Users/gianlucamastrantonio/Politecnico di Torino Staff Dropbox/Gianluca Mastrantonio/lavori/gitrepo/toroidal_projected_normal/functions/general_functions.R")


    # ! new crps

    # new_n_miss <- floor(nrow(theta) * 0.1 * ncol(theta))
    # new_set_miss <- data.frame(obs = 1:nrow(theta), d = rep(1:ncol(theta), each = nrow(theta)))
    # new_index_miss <- sample(1:nrow(new_set_miss), new_n_miss, replace = FALSE)
    # new_index_miss <- new_index_miss[order(new_index_miss)]
    # new_set_miss <- new_set_miss[new_index_miss, ]


    # new_theta_miss <- c()
    # for (ii in 1:nrow(new_set_miss))
    # {
    #  iobs <- new_set_miss[ii, 1]
    #  id <- new_set_miss[ii, 2]
    #  new_theta_miss[ii] <- theta[iobs, id]
    # }


    # nsim_mcmc <- nrow(chain[[1]]$sigma_c_out)
    # for (icol in 1:23)
    # {
    #  for (ichain in 1:length(chain))
    #  {
    #    chain[[ichain]]$crps_tot <- matrix(NA, nrow = nrow(theta_miss), ncol = ncol(theta_miss))
    #    sim <- NA
    #    for (isim in 1:nsim_mcmc)
    #    {
    #      if (ex_df$model == "cwc") {
    #        u <- runif(1, 0, 1)
    #        sim[isim] <- q_wc(u, 0, chain[[ichain]]$prec[isim, icol]) + chain[[ichain]]$circ_mean[isim, icol]
    #      } else {
    #        xs <- rnorm(1, 0, 1)
    #        xc <- rnorm(1, chain[[ichain]]$prec[isim, icol], 1)
    #        sim[isim] <- atan2(xs, xc) + chain[[ichain]]$circ_mean[isim, icol]
    #      }
    #    }
    #    chain[[ichain]]$crps_tot[imiss, icol] <- crps_circ(theta_obs_miss, sim)

    #    # chain[[ichain]]$crps[isim] <-
    #  }
    # }
    # nsim_mcmc <- nrow(chain[[1]]$sigma_c_out)
    # for (imiss in 1:nrow(theta_miss))
    # {
    #  print(imiss)
    #  for (icol in 1:23)
    #  {
    #    theta_obs_miss <- theta_miss[imiss, icol]
    #    for (ichain in 1:length(chain))
    #    {
    #      chain[[ichain]]$crps_tot <- matrix(NA, nrow = nrow(theta_miss), ncol = ncol(theta_miss))
    #      sim <- NA
    #      for (isim in 1:nsim_mcmc)
    #      {
    #        if (ex_df$model == "cwc") {
    #          u <- runif(1, 0, 1)
    #          sim[isim] <- q_wc(u, 0, chain[[ichain]]$prec[isim, icol]) + chain[[ichain]]$circ_mean[isim, icol]
    #        } else {
    #          xs <- rnorm(1, 0, 1)
    #          xc <- rnorm(1, chain[[ichain]]$prec[isim, icol], 1)
    #          sim[isim] <- atan2(xs, xc) + chain[[ichain]]$circ_mean[isim, icol]
    #        }
    #      }
    #      chain[[ichain]]$crps_tot[imiss, icol] <- crps_circ(theta_obs_miss, sim)

    #      # chain[[ichain]]$crps[isim] <-
    #    }
    #  }
    # }


    # if (ex_df$model == cwc) {
    #  for (ichain in 1:length(chain))
    #  {
    #    chain[[ichain]]$crps <- NA
    #    sim <- NA
    #    for (isim in 1:nsim_mcmc)
    #    {
    #      u <- runif(1, 0, 1)
    #      sim[isim] <- q_wc(u, 0, chain[[ichain]]$prec[isim, icol]) + chain[[ichain]]$circ_mean[isim, icol]
    #    }
    #    # chain[[ichain]]$crps[isim] <-
    #  }
    # } else {}


    d <- dim(chain[[1]]$circ_mean)[2]
    list_ret[[id_real]] <- list()
    ## chain 1
    nsim <- nrow(chain[[1]]$sigma_s_out)

    if (ex_df$ind[1] == FALSE) {
      index_non_1 <- which(lower.tri(matrix(0, d, d)))
      # paramaters_cw[[h]] <- list(Sigma_s = res_list$Sigma_s, Sigma_c = res_list$Sigma_c, rho = res_list$rho, mu = res_list$mu, pos_def = res_list$pos_def, index_non_1 = index_non_1)
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

      list_ret[[id_real]]$rhat_sigma <- cbind(rhat_vec[index_non_1], 1)
      list_ret[[id_real]]$ess_sigma <- ess_vec[index_non_1]
      list_ret[[id_real]]$ess_bulk_sigma <- ess_bulk_vec[index_non_1]
      list_ret[[id_real]]$ess_tail_sigma <- ess_tail_vec[index_non_1]

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
      # list_ret[[id_real]]$multiESS_sigma <- multiESS_base2(sim_par[, index_non_1], max_lag = 100)$mESS
      # list_ret[[id_real]]$multiESS_sigma_stand <- list_ret[[id_real]]$multiESS_sigma / ((d^2 - d) / 2)
      # sim_par <- rbind(chain1$sigma_s_out, chain2$sigma_s_out)

      # true_par <- matrix(c(chain[[1]]$Sigma_s), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
      # centered_par <- sim_par - true_par
      ## R x d matrix
      # bias <- colMeans(sim_par - true_par)
      # mse <- colMeans((sim_par - true_par)^2)

      qq1 <- apply(sim_par, 2, quantile, probs = 0.025)
      qq3 <- apply(sim_par, 2, quantile, probs = 0.975)
      # coverage <- (qq1 <= 0) & (qq3 >= 0)

      list_ret[[id_real]]$ICL_sigma <- qq3 - qq1

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
      list_ret[[id_real]]$rhat_sigma_c <- cbind(rhat_vec[index_non_1], 1)
      list_ret[[id_real]]$ess_sigma_c <- ess_vec[index_non_1]
      list_ret[[id_real]]$ess_bulk_sigma_c <- ess_bulk_vec[index_non_1]
      list_ret[[id_real]]$ess_tail_sigma_c <- ess_tail_vec[index_non_1]
      # mcmc_list <- mcmc.list(
      #  mcmc(chain1$sigma_c_out),
      #  mcmc(chain2$sigma_c_out)
      # )
      # rhat <- gelman.diag(mcmc_list, multivariate = FALSE)
      # list_ret[[id_real]]$rhat_sigma_c <- rhat$psrf[index_non_1,]
      # rhsat_sigma_c <- rhat$psrf


      # ess <- effectiveSize(mcmc_list)

      # list_ret[[id_real]]$ess_sigma_c <- ess[index_non_1]

      sim_par <- do.call(
        rbind,
        lapply(chain, function(ch) ch$sigma_c_out)
      )
      # list_ret[[id_real]]$multiESS_sigma_c <- multiESS_base2(sim_par[, index_non_1], max_lag = 100)$mESS
      # list_ret[[id_real]]$multiESS_sigma_c_stand <- list_ret[[id_real]]$multiESS_sigma_c / ((d^2 - d) / 2)

      # true_par <- matrix(c(chain[[1]]$Sigma_c), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
      # centered_par <- sim_par - true_par
      ## R x d matrix
      # bias <- colMeans(sim_par - true_par)
      # mse <- colMeans((sim_par - true_par)^2)

      qq1 <- apply(sim_par, 2, quantile, probs = 0.025)
      qq3 <- apply(sim_par, 2, quantile, probs = 0.975)
      coverage <- (qq1 <= 0) & (qq3 >= 0)

      # list_ret[[id_real]]$bias_sigma_c <- bias[index_non_1]
      # list_ret[[id_real]]$mse_sigma_c <- mse[index_non_1]
      # list_ret[[id_real]]$coverage_sigma_c <- coverage[index_non_1]
      list_ret[[id_real]]$ICL_sigma_c <- qq3 - qq1
    }


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
    chains <- lapply(chain, function(ch) ch[["prec"]])

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

    list_ret[[id_real]]$rhat_prec <- cbind(rhat_vec, 1)
    list_ret[[id_real]]$ess_prec <- ess_vec
    list_ret[[id_real]]$ess_bulk_prec <- ess_bulk_vec
    list_ret[[id_real]]$ess_tail_prec <- ess_tail_vec

    sim_par <- do.call(
      rbind,
      lapply(chain, function(ch) ch$prec)
    )

    # true_par <- matrix(c(chain[[1]]$rho), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
    # centered_par <- sim_par - true_par
    ## R x d matrix
    # bias <- colMeans(sim_par - true_par)
    # mse <- colMeans((sim_par - true_par)^2)

    qq1 <- apply(sim_par, 2, quantile, probs = 0.025)
    qq3 <- apply(sim_par, 2, quantile, probs = 0.975)
    # coverage <- (qq1 <= 0) & (qq3 >= 0)

    # list_ret[[id_real]]$bias_rho <- bias
    # list_ret[[id_real]]$mse_rho <- mse
    # list_ret[[id_real]]$coverage_rho <- coverage
    list_ret[[id_real]]$ICL_prec <- qq3 - qq1

    # ! mu

    # tot_mean <- rbind(chain1$mu_out, chain2$mu_out)
    tot_mean <- do.call(
      rbind,
      lapply(chain, function(ch) ch$circ_mean)
    )
    circ_mean_vec <- c()
    for (j in 1:d)
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
    # list_ret[[id_real]]$rhat_mu <- rhat$psrf


    # ess <- effectiveSize(mcmc_list)

    # list_ret[[id_real]]$ess_mu <- ess
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

    list_ret[[id_real]]$rhat_mu <- cbind(rhat_vec, 1)
    list_ret[[id_real]]$ess_mu <- ess_vec
    list_ret[[id_real]]$ess_bulk_mu <- ess_bulk_vec
    list_ret[[id_real]]$ess_tail_mu <- ess_tail_vec

    sim_par <- do.call(
      rbind,
      lapply(chain, function(ch) ch$mu_app)
    )

    # true_par <- matrix(c((mu - circ_mean_vec + pi) %% (2 * pi)), ncol = ncol(sim_par), nrow = nrow(sim_par), byrow = TRUE)
    ## centered_par <- sim_par - true_par
    # centered_par <- atan2(sin(sim_par - true_par), cos(sim_par - true_par))
    ## R x d matrix
    # bias <- colMeans(centered_par)
    # mse <- colMeans(centered_par^2)

    qq1 <- apply(sim_par, 2, quantile, probs = 0.025)
    qq3 <- apply(sim_par, 2, quantile, probs = 0.975)
    # coverage <- (qq1 <= 0) & (qq3 >= 0)

    # list_ret[[id_real]]$bias_mu <- bias
    # list_ret[[id_real]]$mse_mu <- mse
    # list_ret[[id_real]]$coverage_mu <- coverage
    list_ret[[id_real]]$ICL_mu <- qq3 - qq1

    # !crps
    crps_circ <- function(y, x) {
      cx <- cos(x)
      sx <- sin(x)
      L <- length(x)

      mean_d1 <- 1 - cos(y) * mean(cx) - sin(y) * mean(sx)

      R2 <- sum(cx)^2 + sum(sx)^2
      mean_dx <- 1 - R2 / L^2

      mean_d1 - 0.5 * mean_dx
    }
    na_list <- list()
    theta_no_na <- theta
    for (id in 1:d)
    {
      na_list[[id]] <- which(is.na(theta_no_na[, id]))
    }
    for (ichain in 1:length(chain))
    {
      chain[[ichain]]$crps <- matrix(NA, nrow(index_miss), ncol = 1)
      out_miss <- chain[[ichain]]$missig_out
      hh <- 1
      for (id in 1:d)
      {
        for (im in 1:length(na_list[[id]]))
        {
          chain[[ichain]]$crps[hh] <- crps_circ(theta_all[na_list[[id]][im], id], out_miss[[id]][, im])
          hh <- hh + 1
        }
      }
    }
    chains <- lapply(chain, function(ch) ch[["crps"]])

    # x <- do.call(abind::abind, c(chains, along = 3))
    # x <- aperm(x, c(1, 3, 2))
    # dimnames(x) <- list(
    #  NULL,
    #  paste0("chain", 1:length(chain)),
    #  paste0("par", 1:dim(x)[3])
    # )

    # sim_par <- do.call(
    #  rbind,
    #  lapply(chain, function(ch) ch$crps)
    # )
    list_ret[[id_real]]$crps <- mean(unlist(chains))

    # !crps new

    # chains <- lapply(chain, function(ch) ch[["crps_tot_in"]])

    # x <- do.call(abind::abind, c(chains, along = 3))
    # x <- aperm(x, c(1, 3, 2))
    # dimnames(x) <- list(
    #  NULL,
    #  paste0("chain", 1:length(chain)),
    #  paste0("par", 1:dim(x)[3])
    # )

    # sim_par <- do.call(
    #  rbind,
    #  lapply(chain, function(ch) ch$crps_tot_in)
    # )
    # list_ret[[id_real]]$crps_tot_in <- mean(sim_par)


    # ! PLOTS

    #pdf(paste(dir_out, "Ind=", ex_df$ind[1], "_do_small=", ex_df$do_small[1], "_ntry=", ex_df$ntry[1], "_molt_iter=", ex_df$molt_iter[1], "_model=", ex_df$model[1], "setting_id=", ex_df$setting_id[1], ".pdf", sep = ""), width = 15, height = 15)
    #par(mfrow = c(3, 3))
    #for (id in 1:d)
    #{
    #  plot(chain[[1]]$mu_app[, id], type = "l")
    #  for (iii in 2:length(chain))
    #  {
    #    lines(chain[[iii]]$mu_app[, id], col = iii)
    #  }

    #  plot(density(chain[[1]]$mu_app[, id]))
    #  for (iii in 2:length(chain))
    #  {
    #    lines(density(chain[[iii]]$mu_app[, id]), col = iii)
    #  }
    #  abline(v = pi, col = 2)
    #  plot_multi_acf(chain, "mu_app", id, max_lag = 20)
    #}
    #par(mfrow = c(3, 3))
    #for (id in 1:d) {
    #  plot(chain[[1]]$prec[, id], type = "l")
    #  for (iii in 2:length(chain)) {
    #    lines(chain[[iii]]$prec[, id], col = iii)
    #  }


    #  plot(density(chain[[1]]$prec[, id]))
    #  for (iii in 2:length(chain)) {
    #    lines(density(chain[[iii]]$prec[, id]), col = iii)
    #  }


    #  plot_multi_acf(chain, "prec", id, max_lag = 20)
    #}

    #if (ex_df$ind[1] == FALSE) {
    #  hhhh <- 1
    #  for (id in 1:d) {
    #    for (jd in 1:d) {
    #      if (jd < id) {
    #        plot(chain[[1]]$sigma_s_out[, hhhh],
    #          type = "l"
    #        )
    #        for (iii in 2:length(chain)) {
    #          lines(chain[[iii]]$sigma_s_out[, hhhh], col = iii)
    #        }


    #        plot(density(chain[[1]]$sigma_s_out[, hhhh]),
    #          main = round(rhsat_sigma[hhhh, 1], 3)
    #        )
    #        for (iii in 2:length(chain)) {
    #          lines(density(chain[[iii]]$sigma_s_out[, hhhh]), col = iii)
    #        }


    #        plot_multi_acf(chain, "sigma_s_out", hhhh, max_lag = 20, par_title = paste("Sigma_s[", id, ",", jd, "]"))
    #      }


    #      hhhh <- hhhh + 1
    #    }
    #  }

    #  par(mfrow = c(3, 3))
    #  hhhh <- 1
    #  for (id in 1:d) {
    #    for (jd in 1:d) {
    #      if (jd < id) {
    #        plot(chain[[1]]$sigma_c_out[, hhhh],
    #          type = "l"
    #        )
    #        for (iii in 2:length(chain)) {
    #          lines(chain[[iii]]$sigma_c_out[, hhhh], col = iii)
    #        }


    #        plot(density(chain[[1]]$sigma_c_out[, hhhh]))
    #        for (iii in 2:length(chain)) {
    #          lines(density(chain[[iii]]$sigma_c_out[, hhhh]), col = iii)
    #        }


    #        plot_multi_acf(chain, "sigma_c_out", hhhh, max_lag = 20, par_title = paste("Sigma_c[", id, ",", jd, "]"))
    #      }


    #      hhhh <- hhhh + 1
    #    }
    #  }
    #}
    #par(mfrow = c(3, 3))

    #dev.off()

    save(chain, file = paste(dir_out, "Ind=", ex_df$ind[1], "_do_small=", ex_df$do_small[1], "_ntry=", ex_df$ntry[1], "_molt_iter=", ex_df$molt_iter[1], "_model=", ex_df$model[1], "setting_id=", ex_df$setting_id[1], ".RData", sep = ""))
  }
}

par_ret <- df[set_id_unique, ]

save(par_ret, list_ret, file = "real data/output_diagnostic/ Realnew_post_analisi_results_stat.Rdata")
##
# SECTION:
# data_res_pn_sigma_ind <- matrix(NA, ncol = 3, nrow = 3)
# data_res_pn_sigma_dep <- matrix(NA, ncol = 3, nrow = 3)

# data_res_pn_mu_ind <- matrix(NA, ncol = 3, nrow = 3)
# data_res_pn_mu_dep <- matrix(NA, ncol = 3, nrow = 3)

# data_res_pn_kappa_ind <- matrix(NA, ncol = 3, nrow = 3)
# data_res_pn_kappa_dep <- matrix(NA, ncol = 3, nrow = 3)

# data_an_sp <- data.frame(init_sel = rep(NA, length(list_ret_pn)), d_sel = NA, n_sel = NA, kappa_sel = NA, select_sigma = NA, seed_data = NA, ess_sigma = NA, ess_mu = NA, ess_kappa = NA, rhat_sigma = NA, rhat_mu = NA, rhat_kappa = NA)

# for (ii in 1:length(list_ret_pn))
# {
#  data_an_sp$init_sel[ii] <- list_ret_pn[[ii]]$parameters$init_sel
#  data_an_sp$d_sel[ii] <- list_ret_pn[[ii]]$parameters$d_sel
#  data_an_sp$n_sel[ii] <- list_ret_pn[[ii]]$parameters$n_sel
#  data_an_sp$kappa_sel[ii] <- list_ret_pn[[ii]]$parameters$kappa_sel
#  data_an_sp$select_sigma[ii] <- list_ret_pn[[ii]]$parameters$select_sigma
#  data_an_sp$seed_data[ii] <- list_ret_pn[[ii]]$parameters$seed_data
#  data_an_sp$ess_sigma[ii] <- min(list_ret_pn[[ii]]$ess_sigma)
#  data_an_sp$ess_mu[ii] <- min(list_ret_pn[[ii]]$ess_mu)
#  data_an_sp$ess_kappa[ii] <- min(list_ret_pn[[ii]]$ess_kappa)
#  data_an_sp$rhat_sigma[ii] <- max(list_ret_pn[[ii]]$rhat_sigma)
#  data_an_sp$rhat_mu[ii] <- max(list_ret_pn[[ii]]$rhat_mu)
#  data_an_sp$rhat_kappa[ii] <- max(list_ret_pn[[ii]]$rhat_kappa)
# }


#par_ret$rhat_sigma <- NA
#par_ret$rhat_sigma_c <- NA
#par_ret$rhat_mu <- NA
#par_ret$rhat_prec <- NA
#par_ret$ess_sigma <- NA
#par_ret$ess_sigma_c <- NA
#par_ret$ess_mu <- NA
#par_ret$ess_prec <- NA
#par_ret$crps <- NA
#par_ret$crps_tot_in <- NA

# rhat_sigma = NA, rhat_sigma_c = NA, rhat_mu = NA, rhat_kappa = NA, ess_sigma = NA, ess_sigma_c = NA, ess_mu = NA, ess_kappa = NA, coverage_sigma = NA, coverage_sigma_c = NA, coverage_mu = NA, coverage_kappa = NA

# <- data.frame(init_sel = rep(NA, length(list_ret_cw)), d_sel = NA, n_sel = NA, kappa_sel = NA, select_sigma = NA, seed_data = NA, rhat_sigma = NA, rhat_sigma_c = NA, rhat_mu = NA, rhat_kappa = NA, ess_sigma = NA, ess_sigma_c = NA, ess_mu = NA, ess_kappa = NA, coverage_sigma = NA, coverage_sigma_c = NA, coverage_mu = NA, coverage_kappa = NA)
# for (ii in 1:length(list_ret_cw)) {
#  bias_sigma <- list()
# }
# bias_sigma_c <- list()
# bias_mu <- list()
# bias_prec <- list()

# mse_sigma <- list()
# mse_sigma_c <- list()
# mse_mu <- list()
# mse_prec <- list()

#icl_sigma <- list()
#icl_sigma_c <- list()
#icl_mu <- list()
#icl_prec <- list()


## sima_pos_prop <- rep(NA, length(par_ret))
#for (ii in 1:dim(par_ret)[1])
#{
#  # list_ret_pn[[ii]]$ess_kappa <- list_ret_pn[[ii]]$ess_rho
#  # list_ret_pn[[ii]]$rhat_kappa <- list_ret_pn[[ii]]$rhat_rho

#  if (par_ret$ind[ii] == FALSE) {
#    par_ret$rhat_sigma[ii] <- tex_quantile(1 - c(0.025, 0.001), list_ret[[ii]]$rhat_sigma[, 1])
#    par_ret$rhat_sigma_c[ii] <- tex_quantile(1 - c(0.025, 0.001), list_ret[[ii]]$rhat_sigma_c[, 1])

#    par_ret$ess_sigma[ii] <- tex_quantile(c(0.001, 0.025), list_ret[[ii]]$ess_sigma)
#    par_ret$ess_sigma_c[ii] <- tex_quantile(c(0.001, 0.025), list_ret[[ii]]$ess_sigma_c)
#  }
#  par_ret$rhat_mu[ii] <- tex_quantile(1 - c(0.025, 0.001), list_ret[[ii]]$rhat_mu[, 1])
#  par_ret$rhat_prec[ii] <- tex_quantile(1 - c(0.025, 0.001), list_ret[[ii]]$rhat_prec[, 1])


#  par_ret$ess_mu[ii] <- tex_quantile(c(0.001, 0.025), list_ret[[ii]]$ess_mu)
#  par_ret$ess_prec[ii] <- tex_quantile(c(0.001, 0.025), list_ret[[ii]]$ess_prec)
#  par_ret$crps[ii] <- mean(list_ret[[ii]]$crps)
#  par_ret$crps_tot_in[ii] <- mean(list_ret[[ii]]$crps_tot_in)

#  # sima_pos_prop[ii] <- mean(paramaters_pn[[ii]]$pos_def[[1]])


#  icl_sigma[[ii]] <- list_ret[[ii]]$ICL_sigma
#  icl_sigma_c[[ii]] <- list_ret[[ii]]$ICL_sigma_c
#  icl_mu[[ii]] <- list_ret[[ii]]$ICL_mu
#  # list_ret[[ii]]$ICL_prec <- list_ret[[ii]]$ICL_rho
#  # icl_prec[[ii]] <- list_ret[[ii]]$ICL_prec
#  icl_prec[[ii]] <- list_ret[[ii]]$ICL_prec
#}


#par_ret


#par_ret[par_ret$only_ess == TRUE & par_ret$do_best_init == TRUE & par_ret$ntry == 60, ]


#par_ret[par_ret$ind == T, ]
#par_ret[par_ret$ind == FALSE & par_ret$only_ess == T, ]
#par_ret[par_ret$ind == FALSE & par_ret$only_ess == F, ]

#par_ret[par_ret$ind == FALSE & par_ret$do_small == 10, ]
#par_ret[par_ret$ind == FALSE & par_ret$do_small == 15, ]
#par_ret[par_ret$ind == FALSE & par_ret$do_small == 20, ]
#par_ret[par_ret$ind == FALSE & par_ret$do_small == 25, ]

## SECTION:
#ntry <- 80
#molt_iter <- 12
#wmodel <- "tpn"
#wmodel <- "cwc"

#rownames(Sigma) <- c("\\hat{R}", "ESS", "CIL")
#rownames(Sigma_c) <- c("\\hat{R}", "ESS", "CIL")
#rownames(Mu) <- c("\\hat{R}", "ESS", "CIL")
#rownames(Kappa) <- c("\\hat{R}", "ESS", "CIL")

## w_post <- which(par_ret$ind == ind & par_ret$ntry == ntry & par_ret$do_small == do_small & par_ret$molt_iter == molt_iter & par_ret$model == model)
## w_data <- which(df$ind == ind & df$ntry == ntry & df$do_small == do_small & df$molt_iter == molt_iter & df$model == model)
## name_file <- ff[w_data]
## data_plot <- par_ret[w_post, ]

#w_data <- which(par_ret$only_ess == TRUE & par_ret$do_best_init == TRUE & par_ret$ntry == 60 & par_ret$molt_iter == 2)

#par_ret_plot <- par_ret[w_data, ]

## par_ret[par_ret$only_ess == TRUE & par_ret$do_best_init == TRUE & par_ret$ntry == 60 & par_ret$molt_iter == 2 & par_ret$model == "tpn", ]


#wplot1 <- which(par_ret_plot$model == "tpn" & par_ret_plot$ind == TRUE)
#wplot2 <- which(par_ret_plot$model == "cwc" & par_ret_plot$ind == TRUE)
#wplot3 <- which(par_ret_plot$model == "tpn" & par_ret_plot$ind == FALSE)
#wplot4 <- which(par_ret_plot$model == "cwc" & par_ret_plot$ind == FALSE)

#a_res <- par_ret_plot[wplot1, ]
#b_res <- par_ret_plot[wplot2, ]
#c_res <- par_ret_plot[wplot3, ]
#d_res <- par_ret_plot[wplot4, ]

#Sigma[1, ] <- c("", "", c_res$rhat_sigma, d_res$rhat_sigma)
#Sigma[2, ] <- c("", "", c_res$ess_sigma, d_res$ess_sigma)
#Sigma[3, ] <- c(
#  tex_quantile(c(0.025, 0.975), icl_sigma[[wplot1]]),
#  tex_quantile(c(0.025, 0.975), icl_sigma[[wplot2]]),
#  tex_quantile(c(0.025, 0.975), icl_sigma[[wplot3]]),
#  tex_quantile(c(0.025, 0.975), icl_sigma[[wplot4]])
#)
#Sigma_c[1, ] <- c("", "", c_res$rhat_sigma_c, d_res$rhat_sigma_c)
#Sigma_c[2, ] <- c("", "", c_res$ess_sigma_c, d_res$ess_sigma_c)
#Sigma_c[3, ] <- c(
#  tex_quantile(c(0.025, 0.975), icl_sigma_c[[wplot1]]),
#  tex_quantile(c(0.025, 0.975), icl_sigma_c[[wplot2]]),
#  tex_quantile(c(0.025, 0.975), icl_sigma_c[[wplot3]]),
#  tex_quantile(c(0.025, 0.975), icl_sigma_c[[wplot4]])
#)


## MU
#Mu[1, ] <- c(a_res$rhat_mu, b_res$rhat_mu, c_res$rhat_mu, d_res$rhat_mu)
#Mu[2, ] <- c(a_res$ess_mu, b_res$ess_mu, c_res$ess_mu, d_res$ess_mu)
#Mu[3, ] <- c(
#  tex_quantile(c(0.025, 0.975), icl_mu[[wplot1]]),
#  tex_quantile(c(0.025, 0.975), icl_mu[[wplot2]]),
#  tex_quantile(c(0.025, 0.975), icl_mu[[wplot3]]),
#  tex_quantile(c(0.025, 0.975), icl_mu[[wplot4]])
#)


## KAPPA
#Kappa[1, ] <- c(a_res$rhat_prec, b_res$rhat_prec, c_res$rhat_prec, d_res$rhat_prec)
#Kappa[2, ] <- c(a_res$ess_prec, b_res$ess_prec, c_res$ess_prec, d_res$ess_prec)

#Kappa[3, ] <- c(
#  tex_quantile(c(0.025, 0.975), icl_prec[[wplot1]]),
#  tex_quantile(c(0.025, 0.975), icl_prec[[wplot2]]),
#  tex_quantile(c(0.025, 0.975), icl_prec[[wplot3]]),
#  tex_quantile(c(0.025, 0.975), icl_prec[[wplot4]])
#)

#tt <- table_tex(Sigma, Sigma_c, Mu, Kappa,
#  groups = c("Identity", "Full"),
#  n_vals = c(25, 50, 25, 50),
#  caption = "PN",
#  # prec_name = "\\boldsymbol{\\kappa}"
#  prec_name = "\\boldsymbol{\\lambda}"
#)

#cat(paste(tt, collapse = "\n"))


## SECTION: new CRPS

#### comptue crps
## crps_circ <- function(real_data, missing_vec) {
##  dd <- c(real_data, missing_vec)

##  dist_mat <- 1 - cos(as.matrix(dist(dd)))
##  L <- length(missing_vec)
##  return(sum(dist_mat[1, -1]) / L - 1 / (2 * L^2) * sum(c(dist_mat[-1, -1])))
## }


#w_sel_small <- which(par_ret$ntry == 40 & par_ret$molt_iter == TRUE & par_ret$do_small == TRUE)
#w_sel_large <- which(par_ret$ntry == 40 & par_ret$molt_iter == TRUE & par_ret$do_small == FALSE)
#res_model_small <- par_ret[w_sel_small, ]
#res_model_large <- par_ret[w_sel_large, ]
#best

#chosen_model <- w_sel_small[1]

#w_which_df <- df$setting_id == chosen_model
#ex_df <- df[which(w_which_df)[1], ]

#ff_plot <- ff[w_which_df]
#set_id_unique[id_real] <- c(which(w_which_df)[1])
## print(id_real)
## print(ex_df)

#chain <- list()


#load(paste(dir_data, ff_plot[1], sep = ""))
#if (ex_df$model[1] == "tpn") {
#  prec <- kappa_out
#} else {
#  prec <- rho_out
#}
#res_list <- list(
#  circ_mean = mu_out,
#  prec = prec,
#  sigma_s_out = sigma_s_out,
#  sigma_c_out = sigma_c_out
#)
#chain[[1]] <- res_list
#load(paste(dir_data, ff_plot[2], sep = ""))

#if (ex_df$model[1] == "tpn") {
#  prec <- kappa_out
#} else {
#  prec <- rho_out
#}
#res_list <- list(
#  circ_mean = mu_out,
#  prec = prec,
#  sigma_s_out = sigma_s_out,
#  sigma_c_out = sigma_c_out
#)
#chain[[2]] <- res_list

#if (length(ff_plot) >= 3) {
#  load(paste(dir_data, ff_plot[3], sep = ""))

#  if (ex_df$model[1] == "tpn") {
#    prec <- kappa_out
#  } else {
#    prec <- rho_out
#  }
#  res_list <- list(
#    circ_mean = mu_out,
#    prec = prec,
#    sigma_s_out = sigma_s_out,
#    sigma_c_out = sigma_c_out
#  )
#  chain[[3]] <- res_list
#}
#if (length(ff_plot) >= 4) {
#  load(paste(dir_data, ff_plot[4], sep = ""))

#  if (ex_df$model[1] == "tpn") {
#    prec <- kappa_out
#  } else {
#    prec <- rho_out
#  }
#  res_list <- list(
#    circ_mean = mu_out,
#    prec = prec,
#    sigma_s_out = sigma_s_out,
#    sigma_c_out = sigma_c_out
#  )
#  chain[[4]] <- res_list
#}
#if (length(ff_plot) >= 5) {
#  load(paste(dir_data, ff_plot[5], sep = ""))

#  if (ex_df$model[1] == "tpn") {
#    prec <- kappa_out
#  } else {
#    prec <- rho_out
#  }
#  res_list <- list(
#    circ_mean = mu_out,
#    prec = prec,
#    sigma_s_out = sigma_s_out,
#    sigma_c_out = sigma_c_out
#  )
#  chain[[5]] <- res_list
#}
### parameters
#chains <- lapply(chain, function(ch) ch[["sigma_s_out"]])
#x <- do.call(abind::abind, c(chains, along = 3))
#x <- aperm(x, c(1, 3, 2))
#dimnames(x) <- list(
#  NULL,
#  paste0("chain", 1:length(chain)),
#  paste0("par", 1:dim(x)[3])
#)

#est_sigma_s <- do.call(
#  rbind,
#  lapply(chain, function(ch) ch$sigma_s_out)
#)


#chains <- lapply(chain, function(ch) ch[["circ_mean"]])
#x <- do.call(abind::abind, c(chains, along = 3))
#x <- aperm(x, c(1, 3, 2))
#dimnames(x) <- list(
#  NULL,
#  paste0("chain", 1:length(chain)),
#  paste0("par", 1:dim(x)[3])
#)

#est_mu <- do.call(
#  rbind,
#  lapply(chain, function(ch) ch$circ_mean)
#)


#chains <- lapply(chain, function(ch) ch[["prec"]])
#x <- do.call(abind::abind, c(chains, along = 3))
#x <- aperm(x, c(1, 3, 2))
#dimnames(x) <- list(
#  NULL,
#  paste0("chain", 1:length(chain)),
#  paste0("par", 1:dim(x)[3])
#)

#est_prec <- do.call(
#  rbind,
#  lapply(chain, function(ch) ch$prec)
#)


##
#post_mean_sigma <- matrix(apply(est_sigma_s, 2, mean), nrow = 25)
#post_mean_mu <- apply(est_mu, 2, mean)
#post_mean_prec <- apply(est_prec, 2, mean)
## list_ret[[id_real]]$crps <- sim_par


#coords_info <- metadata_staz_code[w, ]


## plot(coords_info[, 2:3])

#library(ggplot2)
#library(tidyverse)
#library(dplyr)
#library(ggpubr)
#library(gridExtra)
#library(mcclust.ext)

#library(ggplot2)

#library(sf)

#library(rnaturalearth)

#library(rnaturalearthdata)


#library(ggplot2)
#library(sf)
#library(rnaturalearth)
#library(rnaturalearthdata)

## Example: your matrix has longitude in column 1 and latitude in column 2
#coords <- as.data.frame(coords_info[, c(3, 2)])
#names(coords) <- c("lon", "lat")


#coords$mu <- post_mean_mu
#coords$prec <- post_mean_prec

## rescale arrow lengths
#scale_arrow <- 1 / max(coords$prec) * 3

#coords$dx <- cos(coords$mu) * coords$prec * scale_arrow
#coords$dy <- sin(coords$mu) * coords$prec * scale_arrow

#coords$lon_end <- coords$lon + coords$dx
#coords$lat_end <- coords$lat + coords$dy

#n <- nrow(coords)

#edges <- data.frame()

#for (i in 1:(n - 1)) {
#  for (j in (i + 1):n) {
#    edges <- rbind(
#      edges,
#      data.frame(
#        x = coords$lon[i],
#        y = coords$lat[i],
#        xend = coords$lon[j],
#        yend = coords$lat[j],
#        sigma = post_mean_sigma[i, j]
#      )
#    )
#  }
#}
#edges <- subset(edges, sigma > 0.3)
## threshold <- quantile(edges$sigma, 0.9)

## edges <- subset(edges, sigma >= threshold)
#edges$width <- 3 * edges$sigma / max(edges$sigma)


## Convert coordinates to sf object
#points_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

## America map, USA example
#library(ggplot2)
#library(sf)
#library(rnaturalearth)

#america <- ne_countries(
#  continent = c("North America", "South America"),
#  returnclass = "sf"
#)

#lakes <- ne_download(
#  scale = 50,
#  type = "lakes",
#  category = "physical",
#  returnclass = "sf"
#)


#ggplot() +
#  geom_sf(data = america, fill = "gray95", color = "gray20", linewidth = 0.6) +
#  geom_sf(data = lakes, fill = "skyblue", color = "skyblue") +
#  geom_sf(data = points_sf, color = "red", size = 2) +
#  coord_sf(xlim = c(-130, -60), ylim = c(20, 55)) +
#  theme_minimal()

#ggplot() +
#  geom_sf(data = america, fill = "gray95", color = "black", linewidth = 0.6) +
#  geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
#  geom_segment(
#    data = coords,
#    aes(x = lon, y = lat, xend = lon_end, yend = lat_end),
#    arrow = arrow(length = unit(0.15, "cm")),
#    linewidth = 0.5,
#    color = "red"
#  ) +
#  coord_sf(xlim = c(-130, -60), ylim = c(20, 55)) +
#  theme_minimal()


#ggplot() +
#  geom_sf(
#    data = america,
#    fill = "gray95",
#    color = "black",
#    linewidth = 0.6
#  ) +
#  geom_sf(
#    data = lakes,
#    fill = "lightblue",
#    color = "lightblue"
#  ) +
#  geom_segment(
#    data = edges,
#    aes(
#      x = x,
#      y = y,
#      xend = xend,
#      yend = yend,
#      linewidth = width
#    ),
#    color = "steelblue",
#    alpha = 0.6
#  ) +
#  geom_point(
#    data = coords,
#    aes(x = lon, y = lat),
#    color = "red",
#    size = 2
#  ) +
#  scale_linewidth_identity() +
#  coord_sf(xlim = c(-130, -60), ylim = c(20, 55)) +
#  theme_minimal()


#library(ggplot2)
#library(reshape2)
#sigma_df <- reshape2::melt(post_mean_sigma)

#ggplot(
#  sigma_df,
#  aes(x = Var2, y = Var1, fill = value)
#) +
#  geom_tile() +
#  scale_y_reverse() +
#  coord_equal() +
#  scale_fill_gradient2(
#    low = "blue",
#    mid = "white",
#    high = "red",
#    midpoint = 0,
#    limits = c(-1, 1)
#  ) +
#  theme_minimal()


#pdf("/Users/gianlucamastrantonio/Politecnico di Torino Staff Dropbox/Gianluca Mastrantonio/lavori/gitrepo/toroidal_projected_normal/real data/ppp.pdf")
#par(mfrow = c(2, 2))
#for (id in 1:23)
#{
#  hist(c(theta[, id], theta[, id] - 2 * pi, theta[, id] + 2 * pi), main = id, xlim = c(0, 2 * pi), breaks = 36)
#  plot(density(c(theta[, id], theta[, id] - 2 * pi, theta[, id] + 2 * pi), from = 0, to = 2 * pi, adjust = 1 / 3), xlim = c(0, 2 * pi))
#}
#dev.off()

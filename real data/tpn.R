library(CholWishart)

library(ggplot2)
library(dplyr)
library(readr)
library(lubridate)
library(VGAM)
library(MCMCpack)
library(Rfast)
library(MCMCpack)
library(MASS)
library(truncnorm)
library(matrixcalc)
library(LaplacesDemon)
library(toroidalPNcopula)
#source("functions/general_functions.R")
#source("functions/mcmc_tpn.R")
#source("functions/mcmc_tpn_mixture.R")
# source("/beegfs/users/gmastrantonio/tokyo/codes/parameters_mcmc.R")
findmode <- function(x) {
  TT <- table(as.vector(x))
  return(as.numeric(names(TT)[TT == max(TT)][1]))
}
# load("real data/data/data_stations_code.RData")




load("real data/data/gauge.Rdata")



# args <- c(1 ,1 ,1, 1, 1, 80 ,1 ,1, 1)
#args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
do_best_init <- c(T, F)[as.integer(args[2])]
do_small <- c("T1", "T2", "T3", "T4", "T5")[as.integer(args[3])]

do_only_ESS <- c(TRUE, FALSE)[as.integer(args[4])]
type_ess <- as.integer(args[5])
n_test_sigma <- as.integer(args[6])

do_ind <- c(TRUE, FALSE)[as.integer(args[7])]
mmmolt_iter <- c(1, 2)[as.integer(args[8])]

mixt <- c(TRUE, FALSE)[as.integer(args[9])]

Kmax <- as.integer(args[10])

name_sim <- paste(name_sim, "IND", do_ind, "_", do_only_ESS, type_ess, n_test_sigma, "do_best_init=", do_best_init, "do_small=", do_small, "mmmolt_iter=", mmmolt_iter, "mixt=", mixt, "Kmax=", Kmax, sep = "")


# ========
# *SECTION - Functions
# ========
sim_sigma <- function(par1, par2) {
  tryCatch(
    {
      Sigma_s <- riwish(par1, par2)
      chol(abs(Sigma_s))
      return(list(Sigma_s, TRUE))
    },
    error = function(e) {
      return(list(1, FALSE))
    }
  )
}
test_sigma <- function(Sigma_s) {
  tryCatch(
    {
      chol(abs(Sigma_s))
      return(TRUE)
    },
    error = function(e) {
      return(FALSE)
    }
  )
}

crps_circ <- function(real_data, missing_vec) {
  dd <- c(real_data, missing_vec)

  dist_mat <- 1 - cos(as.matrix(dist(dd)))
  L <- length(missing_vec)
  return(sum(dist_mat[1, -1]) / L - 1 / (2 * L^2) * sum(c(dist_mat[-1, -1])))
}

# ========
# * SECTION - data
# ========
# theta <- as.matrix(theta/ 360 * 2 * pi)
# theta_all <- as.matrix(theta_all/ 360 * 2 * pi)
# set.seed(1)
## w <- sample(1:ncol(theta),50)
# if(do_small == "T1")
# {
#  theta <- theta
# }else{
#  error("Invalid value for do_small")
# }


app <- theta
for (id in 1:ncol(theta))
{
  app[which(is.na(app[, id])), id] <- mean(app[, id], na.rm = T)
}


n <- nrow(theta)
d <- ncol(theta)

# * na
set.seed(123)
na_list <- list()
theta_no_na <- theta
for (id in 1:d)
{
  na_list[[id]] <- which(is.na(theta_no_na[, id]))
  if (length(na_list[[id]]) > 0) {
    theta_no_na[na_list[[id]], id] <- runif(length(na_list[[id]]), 0, 2 * pi)
  }
}
theta_no_na <- theta_all
n <- nrow(theta_no_na)
d <- ncol(theta_no_na)
# ========
# * for crps
# ========

# n_miss <- floor(n * 0.1)
# n_miss <- 20
# y_miss <- matrix(NA, nrow = n_miss, ncol = d)
# for (id in 1:d)
# {
#  index_miss <- sample(which(!is.na(theta[, id])), n_miss)
#  y_miss[, id] <- theta[index_miss, id]
#  na_list[[id]] <- c(index_miss, na_list[[id]])
# }


# ========
# * SECTION - MCMC
# ========


mmm <- m_mcmc
set.seed(seed)
if (do_best_init == TRUE) {
  mean_init <- mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x), na.rm = TRUE), sum(cos(x), na.rm = TRUE)))
  kappa_init <- rep(1, d)
  r_init <- matrix(1, n, d)
  x_init <- matrix(0, n, d)
  y_init <- matrix(0, n, d)

  for (i in 1:d)
  {
    theta_app <- theta[, i] - mean_init[i]
    if (sum(is.na(theta_app)) > 0) {
      theta_app[is.na(theta_app)] <- runif(sum(is.na(theta_app)), 0, 2 * pi)
    }

    u <- runif(n)
    C <- mean(cos(theta_app), na.rm = TRUE)
    S <- mean(sin(theta_app), na.rm = TRUE)
    R <- sqrt(C^2 + S^2)
    V <- 1 - R
    kappa_init[i] <- 1 / sqrt(2 * V)


    r_app <- r_rice(n, kappa_init, sigma = 1)


    # r_direct <- sqrt(-2 * log(u))
    # r_app <- r_direct
    # print("AA")
    # print(dim(theta))
    # print(length(r_app))
    # print(length(theta_app))
    x_app <- r_app * cos(theta_app)
    y_app <- r_app * sin(theta_app)

    # x_app <- x_app - mean(x_app) + 1
    # y_app <- y_app - mean(y_app)

    # ttt <- atan2(y_app, x_app)
    # qq <- quantile((ttt) + pi, prob = c(0.15, 0.85)) - pi
    # kappa_init[i] <- cos(qq[2]) / sin(qq[2])
    # x_app <- x_app + kappa_init[i]
    # r_init[, i] <- sqrt(x_app^2 + y_app^2)
    # x_init[, i] <- x_app
    y_init[, i] <- y_app
  }
  if (do_ind == TRUE) {
    sigma_init <- diag(1, d)
  } else {
    sigma_init <- cov(y_init)
    while (inherits(try(chol(abs(sigma_init)), silent = TRUE), "try-error")) {
      sigma_init <- sigma_init + diag(0.01, d)
    }
  }
} else {
  mean_init <- mu_init <- apply(theta, 2, function(x) atan2(sum(sin(x), na.rm = TRUE), sum(cos(x), na.rm = TRUE)))
  kappa_init <- runif(d, 0.7, 0.9)
  r_init <- matrix(runif(n * d, 0.8, 1.2), n, d)
  sigma_init <- diag(1, d)
}
#print("initialization done")
#print(warnings())
start <- Sys.time()

if (mixt == FALSE) {
  out_mcmc <- mcmc_tpn(
    theta = theta_no_na, # the circualr data
    burnin = burnin_mcmc * mmm * mmmolt_iter, # burnin
    thin = thin_mcmc * mmm * mmmolt_iter, # thin
    iterations = iter_mcmc * mmm * mmmolt_iter, # total interations
    # burnin = 10, # burnin
    # thin = 1  , # thin
    # iterations = 30, # total interations
    prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
    prior_mu_var = rep(100000, d),
    prior_kappa_mean = matrix(0, nrow = d, ncol = 1), # the prior for k is  TN(prior_kappa_mean,prior_kappa_var )
    prior_kappa_var = rep(100000, d),
    prior_sigma_nu = d + nu_app, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
    prior_sigma_psi = diag(1, d) * (d + nu_app - d - 1),


    # this section set the initial values of the parameters
    mu_init = mu_init,
    kappa_init = kappa_init,
    sigma_init = sigma_init,
    r_init = r_init,

    # parameters for the adaptive part of Metropolis
    adapt_batch = batch_mcmc,
    adapt_a = a_mcmc,
    adapt_b = b_mcmc,
    adapt_alpha_target = alpha_target,
    sd_mu_scal = 1,
    par_sigma_adapt = par_adapt_mcmc,
    na_index = na_list,
    n_test_sigma = n_test_sigma,
    do_only_ESS = do_only_ESS,
    type_ess = type_ess,
    do_ind = do_ind
  )
  end <- Sys.time()
  runtime <- end - start

  # ========
  # * SECTION - Output
  # ========
  mu_out <- out_mcmc$mu_out
  kappa_out <- out_mcmc$kappa_out
  sigma_s_out <- out_mcmc$sigma_s_out
  sigma_c_out <- out_mcmc$sigma_c_out
  r_out <- out_mcmc$r_out


  # missig_out <- out_mcmc$missig_out
  waic <- out_mcmc$waic
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


  # crps_val <- matrix(0, nrow = n_miss, ncol = d)
  # for (id in 1:d)
  # {
  #  for (imiss in 1:n_miss)
  #  {
  #    crps_val[imiss, id] <- crps_circ(y_miss[imiss, id], missig_out[[id]][, imiss])
  #  }
  # }


  save.image(paste("real data/output/", name_sim, "tpn_seed", seed, ".Rdata", sep = ""))


  #pdf(paste("real data/output/", name_sim, "tpn_seed", seed, ".pdf", sep = ""))


  ## plot(c(crps_val), main = paste( " - ", round(mean(c(waic)), 5)))


  #par(mfrow = c(1, 1))


  #data_plot <- data.frame(var = colMeans(sigma_s_out[, ]), x = rep((1:d), each = d), y = rep((1:d), times = d))
  #p1 <- data_plot %>% ggplot(aes(x = x, y = y, fill = var)) +
  #  geom_tile() +
  #  scale_y_reverse() +
  #  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))
  #print(p1)


  ## for (id in 1:d)
  ## {
  ##  data_plot <- data.frame(var = c(kappa_out[, id, ]), iter = rep(1:dim(kappa_out)[1], times = K), kk = factor(rep(1:K, each = dim(kappa_out)[1])))

  ##  p1 <- data_plot %>% ggplot(aes(x = iter, y = var, col = kk, group = kk)) +
  ##    geom_line() +
  ##    ggtitle(paste("kappa", id))
  ##  print(p1)
  ## }


  #par(mfrow = c(3, 3))
  #for (id in 1:d)
  #{
  #  plot(mu_out[, id], type = "l")
  #}
  #for (id in 1:d)
  #{
  #  plot(kappa_out[, id], type = "l")
  #}

  #h <- 1
  #for (id in 1:d)
  #{
  #  for (jd in 1:d)
  #  {
  #    plot(sigma_s_out[, h], type = "l")
  #    h <- h + 1
  #  }
  #}
  #dev.off()
} else {
  out_mcmc <- mcmc_tpn_mixture(
    theta = theta_no_na, # the circualr data
    burnin = burnin_mcmc * mmm * mmmolt_iter, # burnin
    thin = thin_mcmc * mmm * mmmolt_iter, # thin
    iterations = iter_mcmc * mmm * mmmolt_iter, # total interations
    # burnin = 10, # burnin
    # thin = 1  , # thin
    # iterations = 30, # total interations
    prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
    prior_mu_var = rep(100000, d),
    prior_kappa_mean = matrix(0, nrow = d, ncol = 1), # the prior for k is  TN(prior_kappa_mean,prior_kappa_var )
    prior_kappa_var = rep(100000, d),
    prior_sigma_nu = d + nu_app, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
    prior_sigma_psi = diag(1, d) * (d + nu_app - d - 1),


    # this section set the initial values of the parameters
    mu_init = mu_init,
    kappa_init = kappa_init,
    sigma_init = sigma_init,
    r_init = r_init,

    # parameters for the adaptive part of Metropolis
    adapt_batch = batch_mcmc,
    adapt_a = a_mcmc,
    adapt_b = b_mcmc,
    adapt_alpha_target = alpha_target,
    sd_mu_scal = 1,
    par_sigma_adapt = par_adapt_mcmc,
    na_index = na_list,
    n_test_sigma = n_test_sigma,
    do_only_ESS = do_only_ESS,
    type_ess = type_ess,
    do_ind = do_ind,
    Kmax = Kmax
  )
  end <- Sys.time()
  runtime <- end - start

  # ========
  # * SECTION - Output
  # ========
  mu_out <- out_mcmc$mu_out
  kappa_out <- out_mcmc$kappa_out
  sigma_s_out <- out_mcmc$sigma_s_out
  sigma_c_out <- out_mcmc$sigma_c_out
  r_out <- out_mcmc$r_out
  z_out <- out_mcmc$z_out

  # missig_out <- out_mcmc$missig_out
  waic <- out_mcmc$waic
  ### identification
  mu_out <- mu_out %% (2 * pi)
  nsim <- nrow(mu_out)

  # # # # # # # # # # # # # #
  # The parameters must be indentified
  # # # # # # # # # # # # # #
  for (k in 1:dim(sigma_s_out)[3])
  {
    for (isim in 1:nsim)
    {
      ss <- matrix(sigma_s_out[isim, , k], nrow = d)
      B <- diag(1 / diag(ss)^0.5)

      sigma_s_out[isim, , k] <- B %*% matrix(sigma_s_out[isim, , k], nrow = d) %*% B
      sigma_c_out[isim, , k] <- B %*% matrix(sigma_c_out[isim, , k], nrow = d) %*% B
      kappa_out[isim, , k] <- kappa_out[isim, , k] * diag(B)

      # for (iobs in 1:n)
      # {
      #  r_out[isim, iobs, ] <- r_out[isim, iobs, ] * diag(B)
      # }
    }
  }


  # crps_val <- matrix(0, nrow = n_miss, ncol = d)
  # for (id in 1:d)
  # {
  #  for (imiss in 1:n_miss)
  #  {
  #    crps_val[imiss, id] <- crps_circ(y_miss[imiss, id], missig_out[[id]][, imiss])
  #  }
  # }


  save.image(paste("real data/output/", name_sim, "tpn_seed", seed, ".Rdata", sep = ""))


  #pdf(paste("real data/output/", name_sim, "tpn_seed", seed, ".pdf", sep = ""))


  ## plot(c(crps_val), main = paste( " - ", round(mean(c(waic)), 5)))


  ## pdf(paste("real data/output/",name_sim,"cwc_seed",seed,".pdf", sep = ""))


  #K <- dim(sigma_s_out)[3]
  ## plot(c(crps_val), main = paste(round(mean(c(crps_val)), 5), " - ", round(mean(c(waic)), 5)))

  #k_samp <- rep(NA, dim(z_out)[1])
  #for (isim in 1:dim(z_out)[1])
  #{
  #  k_samp[isim] <- length(unique(z_out[isim, ]))
  #}
  #par(mfrow = c(1, 2))
  #plot(k_samp, type = "l")
  #barplot(table(k_samp))

  #par(mfrow = c(1, 1))
  #zeta_map <- apply(z_out, 2, findmode)
  #par(mfrow = c(1, 1))
  #plot(zeta_map, type = "l")
  #barplot(table(zeta_map))


  #for (k in 1:K)
  #{
  #  data_plot <- data.frame(var = colMeans(sigma_s_out[, , k]), x = rep((1:d), each = d), y = rep((1:d), times = d))
  #  p1 <- data_plot %>% ggplot(aes(x = x, y = y, fill = var)) +
  #    geom_tile() +
  #    scale_y_reverse() +
  #    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))
  #  print(p1)
  #}
  #for (id in 1:d)
  #{
  #  data_plot <- data.frame(var = c(mu_out[, id, ]), iter = rep(1:dim(mu_out)[1], times = K), kk = factor(rep(1:K, each = dim(mu_out)[1])))

  #  p1 <- data_plot %>% ggplot(aes(x = iter, y = var, col = kk, group = kk)) +
  #    geom_line() +
  #    ylim(0, 2 * pi) +
  #    ggtitle(paste("mu", id))
  #  print(p1)
  #}


  #for (id in 1:d)
  #{
  #  data_plot <- data.frame(var = c(kappa_out[, id, ]), iter = rep(1:dim(kappa_out)[1], times = K), kk = factor(rep(1:K, each = dim(kappa_out)[1])))

  #  p1 <- data_plot %>% ggplot(aes(x = iter, y = var, col = kk, group = kk)) +
  #    geom_line() +
  #    ggtitle(paste("kappa", id))
  #  print(p1)
  #}


  #par(mfrow = c(3, 3))
  #for (id in 1:d)
  #{
  #  for (k in 1:K)
  #  {
  #    plot(mu_out[, id, k], type = "l")
  #  }
  #}
  #for (id in 1:d)
  #{
  #  for (k in 1:K)
  #  {
  #    plot(kappa_out[, id, k], type = "l")
  #  }
  #}

  #for (k in 1:K)
  #{
  #  h <- 1
  #  for (id in 1:d)
  #  {
  #    for (jd in 1:d)
  #    {
  #      plot(sigma_s_out[, h, k], type = "l")
  #      h <- h + 1
  #    }
  #  }
  #}
  #dev.off()
}

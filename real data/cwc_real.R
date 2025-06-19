rm(list = ls())
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
source("functions/general_functions.R")
source("functions/mcmc_cwn.R")


load("real data/data/data wind.Rdata")



# ========
# * SECTION - data
# ========
theta <- as.matrix(data_subset_model[, -1] / 350 * 2 * pi)
theta <- theta[, c(1, 2, 3, 5, 4, 6)]

n <- nrow(theta)
d <- ncol(theta)
theta[1, ] <- NA
# * na
na_list <- list()
theta_no_na <- theta
for (id in 1:d)
{
  na_list[[id]] <- which(is.na(theta_no_na[, id]))
  if (length(na_list[[id]]) > 0) {
    theta_no_na[na_list[[id]], id] <- runif(length(na_list[[id]]), 0, 2 * pi)
  }
}

## ========
## * for crps
## ========
set.seed(123)
n_miss <- floor(n * 0.1)
n_miss <- 20
y_miss <- matrix(NA, nrow = n_miss, ncol = d)
for (id in 1:d)
{
  index_miss <- sample(which(!is.na(theta[, id])), n_miss)
  y_miss[, id] <- theta[index_miss, id]
  na_list[[id]] <- c(index_miss, na_list[[id]])
}


# ========
# * SECTION - MCMC
# ========

# sigma_init <- array(NA, c(d, d, K))
# for (k in 1:K)
# {
#  sigma_init[, , k] <- diag(runif(d, 0.5, 1.5))
# }

set.seed(2)
mmm <- 10
out_mcmc <- mcmc_cwc(
  theta = theta_no_na, # the circualr data
  burnin = 1000 * mmm, # burnin
  thin = 1 * mmm, # thin
  iterations = 3000 * mmm, # total interations
  # burnin = 10, # burnin
  # thin = 1, # thin
  # iterations = 20, # total interations
  prior_mu_mean = matrix(0, nrow = d, ncol = 1), # the prior on the mean is N(prior_mu_mean,prior_mu_var )
  prior_mu_var = rep(100000, d),
  prior_rho_a = rep(1, d), # the prior for B()
  prior_rho_b = rep(1, d),
  prior_sigma_nu = d + 2, # the prior for sigma is  IW(prior_sigma_nu, prior_sigma_psi)
  prior_sigma_psi = diag(1, d),


  # this section set the initial values of the parameters
  mu_init = rnorm(d) %% (2 * pi),
  rho_init = abs(runif(d, 0.2, 0.7)),
  sigma_init = diag(1, d),
  r_init = matrix(1, nrow = n, ncol = d),

  # parameters for the adaptive part of Metropolis
  adapt_batch = 50,
  adapt_a = 1000,
  adapt_b = 1200,
  adapt_alpha_target = 0.234,
  sd_mu_scal = 1,
  sd_rho_scal = 0.1,
  par_sigma_adapt = 2000,
  na_index = na_list
)

# ========
# * SECTION - Output
# ========
mu_out <- out_mcmc$mu_out
rho_out <- out_mcmc$rho_out
sigma_s_out <- out_mcmc$sigma_s_out
sigma_c_out <- out_mcmc$sigma_c_out
r_out <- out_mcmc$r_out

missig_out <- out_mcmc$missig_out
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
}


crps_val <- matrix(0, nrow = n_miss, ncol = d)
for (id in 1:d)
{
  for (imiss in 1:n_miss)
  {
    crps_val[imiss, id] <- crps_circ(y_miss[imiss, id], missig_out[[id]][, imiss])
  }
}




save.image(paste("real data/output/cwc.Rdata", sep = ""))



pdf(paste("real data/output/cwc.pdf", sep = ""))



plot(c(crps_val), main = paste(round(mean(c(crps_val)), 5), " - ", round(mean(c(waic)), 5)))



data_plot <- data.frame(var = colMeans(sigma_s_out[, ]), x = rep((1:d), each = d), y = rep((1:d), times = d))
p1 <- data_plot %>% ggplot(aes(x = x, y = y, fill = var)) +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))
print(p1)
# for (id in 1:d)
# {
#  data_plot <- data.frame(var = c(mu_out[, id, ]), iter = rep(1:dim(mu_out)[1], times = K), kk = factor(rep(1:K, each = dim(mu_out)[1])))

#  p1 <- data_plot %>% ggplot(aes(x = iter, y = var, col = kk, group = kk)) +
#    geom_line() +
#    ylim(0, 2 * pi) +
#    ggtitle(paste("mu", id))
#  print(p1)
# }


# for (id in 1:d)
# {
#  data_plot <- data.frame(var = c(rho_out[, id, ]), iter = rep(1:dim(rho_out)[1], times = K), kk = factor(rep(1:K, each = dim(rho_out)[1])))

#  p1 <- data_plot %>% ggplot(aes(x = iter, y = var, col = kk, group = kk)) +
#    geom_line() +
#    ggtitle(paste("rho", id))
#  print(p1)
# }




par(mfrow = c(3, 3))
for (id in 1:d)
{
  plot(mu_out[, id], type = "l")
}
for (id in 1:d)
{
  plot(rho_out[, id], type = "l")
}
h <- 1
for (id in 1:d)
{
  for (jd in 1:d)
  {
    plot(sigma_s_out[, h], type = "l")
    h <- h + 1
  }
}
dev.off()

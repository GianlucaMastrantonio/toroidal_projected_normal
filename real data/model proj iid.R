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
name <- "PN"

setwd("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/Projected Normal on Torus/realdata/")
load("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/Projected Normal on Torus/realdata/dataset/wind/data wind.Rdata")

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
theta <- as.matrix(data_subset_model[, -1] / 350 * 2 * pi)
theta <- theta[, c(1, 2, 3, 5, 4, 6)]


app <- theta
for (id in 1:ncol(theta))
{
  app[which(is.na(app[, id])), id] <- mean(app[, id], na.rm = T)
}


n <- nrow(theta)
d <- ncol(theta)

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

# ========
# * for crps
# ========
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


mmm <- 10
set.seed(1)
# source("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/Projected Normal on Torus/realdata/mcmc_mixture_tpn.R")
source("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/Projected Normal on Torus/simulations/mcmc.R")
out_mcmc <- mcmc_tpn(
  theta = theta_no_na, # the circualr data
  burnin = 1000 * mmm, # burnin
  thin = 1 * mmm, # thin
  iterations = 3000 * mmm, # total interations
  # burnin = 10, # burnin
  # thin = 1, # thin
  # iterations = 20, # total interations
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
  par_sigma_adapt = 2000,
  na_index = na_list
)


# ========
# * SECTION - Output
# ========
mu_out <- out_mcmc$mu_out
kappa_out <- out_mcmc$kappa_out
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
  kappa_out[isim, ] <- kappa_out[isim, ] * diag(B)

  for (iobs in 1:n)
  {
    r_out[isim, iobs, ] <- r_out[isim, iobs, ] * diag(B)
  }
}


crps_val <- matrix(0, nrow = n_miss, ncol = d)
for (id in 1:d)
{
  for (imiss in 1:n_miss)
  {
    crps_val[imiss, id] <- crps_circ(y_miss[imiss, id], missig_out[[id]][, imiss])
  }
}



save.image(paste("plot/tpn.Rdata", sep = ""))



pdf(paste("plot/tpn.pdf", sep = ""))


plot(c(crps_val), main = paste(round(mean(c(crps_val)), 5), " - ", round(mean(c(waic)), 5)))


par(mfrow = c(1, 1))






data_plot <- data.frame(var = colMeans(sigma_s_out[, ]), x = rep((1:d), each = d), y = rep((1:d), times = d))
p1 <- data_plot %>% ggplot(aes(x = x, y = y, fill = var)) +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))
print(p1)



# for (id in 1:d)
# {
#  data_plot <- data.frame(var = c(kappa_out[, id, ]), iter = rep(1:dim(kappa_out)[1], times = K), kk = factor(rep(1:K, each = dim(kappa_out)[1])))

#  p1 <- data_plot %>% ggplot(aes(x = iter, y = var, col = kk, group = kk)) +
#    geom_line() +
#    ggtitle(paste("kappa", id))
#  print(p1)
# }




par(mfrow = c(3, 3))
for (id in 1:d)
{
  plot(mu_out[, id], type = "l")
}
for (id in 1:d)
{
  plot(kappa_out[, id], type = "l")
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

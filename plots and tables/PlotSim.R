library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(gridExtra)

dir_plot <- "plots and tables/out/"

# ========
# PN
# ========
load(paste("simulations/output/tpn_simulations_results -  select_seed= 2 .Rdata", sep = ""))
# load(paste(dir_data, "pn_simulations_parameters -  select_seed= 1 .Rdata", sep = ""))

n_vec <- c(20, 50, 100)
d_vec <- c(6)
w_n <- c()
w_d <- c()
w_k <- c()
for (isim in 1:length(out))
{
  if (out[[isim]]$d %in% d_vec) w_d <- c(w_d, isim)
  if (out[[isim]]$n %in% n_vec) w_n <- c(w_n, isim)
  if (out[[isim]]$kappa[1] == out[[isim]]$kappa[2]) w_k <- c(w_k, isim)
}
w_data <- intersect(w_n, w_d) %>% intersect(w_k)


## create dataset
app <- matrix(seq(1:(6^2)), ncol = 6)
index_cor <- app[lower.tri(app)]


h <- 1
isim <- w_data[1]

out_w <- out[[isim]]

kappa_type <- ifelse(out_w$kappa[3] == 0.49, "1", ifelse(out_w$kappa[3] == 1.1, "2", ifelse(out_w$kappa[3] == out_w$kappa[2], "3", "4")))

nsim_mcmc <- nrow(out_w$mcmc_kappa_out)
data_plot <- data.frame(
  d = out_w$d,
  n = out_w$n,
  # kappa
  kappa_type = kappa_type,
  kappa_sim = matrix(out_w$kappa, ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
  kappa = out_w$kappa_out,
  # sigma
  cor_type = ifelse(out_w$Sigma_s[1, 2] == 0, "1", "2"),
  cor_sim = matrix(out_w$Sigma_s[lower.tri(app)], ncol = length(index_cor), nrow = nsim_mcmc, byrow = T),
  cor = out_w$sigma_s_out[, lower.tri(app)],
  # mu
  mu_type = "1",
  mu_sim = matrix(c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6)[1:out_w$d], ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
  mu = out_w$mu_out,
  index = h
)
for (isim in w_data)
{
  h <- h + 1
  out_w <- out[[isim]]

  kappa_type <- ifelse(out_w$kappa[3] == 0.49, "1", ifelse(out_w$kappa[3] == 1.1, "2", ifelse(out_w$kappa[3] == out_w$kappa[2], "3", "4")))


  nsim_mcmc <- nrow(out_w$mcmc_kappa_out)
  app_frame <- data.frame(
    d = out_w$d,
    n = out_w$n,
    # kappa
    kappa_type = kappa_type,
    kappa_sim = matrix(out_w$kappa, ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
    kappa = out_w$kappa_out,
    # sigma
    cor_type = ifelse(out_w$Sigma_s[1, 2] == 0, "1", "2"),
    cor_sim = matrix(out_w$Sigma_s[lower.tri(app)], ncol = length(index_cor), nrow = nsim_mcmc, byrow = T),
    cor = out_w$sigma_s_out[, lower.tri(app)],
    # mu
    mu_type = "1",
    mu_sim = matrix(c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6)[1:out_w$d], ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
    mu = out_w$mu_out,
    index = h
  )
  data_plot <- rbind(data_plot, app_frame)
}


# ========
# PLOTS
# ========
gg_theme <- theme(
  plot.title = element_text(size = 18), # Increase title size
  axis.title = element_text(size = 16), # Increase axis title size
  axis.text = element_text(size = 20), # Increase axis text size
  legend.title = element_text(size = 25), # Increase legend title size
  legend.text = element_text(size = 18), # Increase legend text size
  legend.position = "bottom",
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  strip.text = element_text(size = 20)
)

for (isigma in 1:2)
{
  #  isigma <- 1
  ## mu



  w_post <- which(colnames(data_plot) %in% paste("mu.", 1:out_w$d, sep = ""))
  w_sim <- which(colnames(data_plot) %in% paste("mu_sim.", 1:out_w$d, sep = ""))

  diff_var <- data_plot[, w_post] - data_plot[, w_sim]
  colnames(diff_var) <- paste("par", 1:out_w$d, sep = "")


  data_plot_app_1 <- data_plot %>%
    mutate(
      var = diff_var
    ) %>%
    unnest(
      cols = var,
      data = .
    ) %>%
    filter(.$cor_type == paste(isigma))

  data_plot_app <- data_plot_app_1 %>%
    select("n", "kappa_type", "cor_type", c(paste("par", 1:out_w$d, sep = ""))) %>%
    pivot_longer(
      cols = paste("par", 1:out_w$d, sep = ""),
      names_to = "type",
      values_to = "samp"
    )

  summary_data <- data_plot_app %>%
    group_by(type, cor_type, n, kappa_type) %>%
    summarize(
      mean = mean(ifelse(samp > pi, samp - 2 * pi, samp)),
      lower = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.025),
      upper = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.975)
    ) %>%
    mutate(kappa = case_when(
      kappa_type == "1" ~ "0.49",
      kappa_type == "2" ~ "1.1",
      kappa_type == "3" ~ "2.25"
    ))


  # summary_data <- summary_data %>%
  #   %>%
  #  mutate(kappa_order = as.numeric(kappa) + as.numeric(kappa_type) / 10)


  pdf(paste(dir_plot, "pn_mu_pos_sigma", c("ind", "dep")[isigma], ".pdf", sep = ""), height = 7, width = 14)
  p <- summary_data %>% ggplot(aes(x = type, y = mean, type = factor(kappa))) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
    geom_errorbar(aes(ymin = lower, ymax = upper, col = kappa), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
    labs(x = "Parameter", y = "95% Credible Interval") +
    gg_theme +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "black", size = 1
    ) +
    scale_x_discrete(labels = c(
      "par1" = expression(mu[1]),
      "par2" = expression(mu[2]),
      "par3" = expression(mu[3]),
      "par4" = expression(mu[4]),
      "par5" = expression(mu[5]),
      "par6" = expression(mu[6])
    )) +
    labs(color = expression(kappa[j])) +
    facet_wrap(~n)
  print(p)
  dev.off()



  ## kappa
  w_post <- which(colnames(data_plot) %in% paste("kappa.", 1:out_w$d, sep = ""))
  w_sim <- which(colnames(data_plot) %in% paste("kappa_sim.", 1:out_w$d, sep = ""))

  diff_var <- data_plot[, w_post] - data_plot[, w_sim]
  colnames(diff_var) <- paste("par", 1:out_w$d, sep = "")
  data_plot_app_1 <- data_plot %>%
    mutate(
      var = diff_var
    ) %>%
    unnest(
      cols = var,
      data = .
    ) %>%
    filter(.$cor_type == paste(isigma))

  data_plot_app <- data_plot_app_1 %>%
    select("n", "kappa_type", "cor_type", c(paste("par", 1:out_w$d, sep = ""))) %>%
    pivot_longer(
      cols = paste("par", 1:out_w$d, sep = ""),
      names_to = "type",
      values_to = "samp"
    )


  # data_plot_app <- cbind(data_plot_app_1, data_plot_app_2) %>% mutate(type = str_replace(type, "mu.", ""))


  summary_data <- data_plot_app %>%
    group_by(type, cor_type, n, kappa_type) %>%
    summarize(
      mean = mean(ifelse(samp > pi, samp - 2 * pi, samp)),
      lower = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.025),
      upper = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.975)
    ) %>%
    mutate(kappa = case_when(
      kappa_type == "1" ~ "0.49",
      kappa_type == "2" ~ "1.1",
      kappa_type == "3" ~ "2.25"
    ))


  pdf(paste(dir_plot, "pn_kappa_pos_sigma", c("ind", "dep")[isigma], ".pdf", sep = ""), height = 7, width = 14)
  p <- summary_data %>% ggplot(aes(x = type, y = mean, type = factor(kappa))) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
    geom_errorbar(aes(ymin = lower, ymax = upper, col = kappa), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
    labs(x = "Parameter", y = "95% Credible Interval") +
    gg_theme +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "black", size = 1
    ) +
    scale_x_discrete(labels = c(
      "par1" = expression(kappa[1]),
      "par2" = expression(kappa[2]),
      "par3" = expression(kappa[3]),
      "par4" = expression(kappa[4]),
      "par5" = expression(kappa[5]),
      "par6" = expression(kappa[6])
    )) +
    labs(color = expression(kappa[j])) +
    facet_wrap(~n)
  print(p)
  dev.off()

  ## rho
  w_post <- which(colnames(data_plot) %in% paste("cor.", 1:length(index_cor), sep = ""))
  w_sim <- which(colnames(data_plot) %in% paste("cor_sim.", 1:length(index_cor), sep = ""))

  diff_var <- data_plot[, w_post] - data_plot[, w_sim]
  colnames(diff_var) <- paste("par", 1:length(index_cor), sep = "")
  data_plot_app_1 <- data_plot %>%
    mutate(
      var = diff_var
    ) %>%
    unnest(
      cols = var,
      data = .
    ) %>%
    filter(.$cor_type == paste(isigma))

  data_plot_app <- data_plot_app_1 %>%
    select("n", "kappa_type", "cor_type", c(paste("par", 1:length(index_cor), sep = ""))) %>%
    pivot_longer(
      cols = paste("par", 1:length(index_cor), sep = ""),
      names_to = "type",
      values_to = "samp"
    )


  # data_plot_app <- cbind(data_plot_app_1, data_plot_app_2) %>% mutate(type = str_replace(type, "mu.", ""))


  summary_data <- data_plot_app %>%
    group_by(type, cor_type, n, kappa_type) %>%
    summarize(
      mean = mean(ifelse(samp > pi, samp - 2 * pi, samp)),
      lower = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.025),
      upper = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.975)
    ) %>%
    mutate(kappa = case_when(
      kappa_type == "1" ~ "0.49",
      kappa_type == "2" ~ "1.1",
      kappa_type == "3" ~ "2.25"
    ))

  # summary_data <- summary_data %>%
  #  mutate(kappa = case_when(
  #    kappa_type == "1" ~ "1",
  #    kappa_type == "2" ~ "10",
  #    kappa_type == "3" ~ "20",
  #    kappa_type == "4" & type == "par1" ~ "1",
  #    kappa_type == "4" & type == "par2" ~ "10",
  #    kappa_type == "4" & type == "par3" ~ "20",
  #    kappa_type == "4" & type == "par4" ~ "1",
  #    kappa_type == "4" & type == "par5" ~ "10",
  #    kappa_type == "4" & type == "par6" ~ "20"
  #  )) %>%
  #  mutate(kappa_order = as.numeric(kappa) + as.numeric(kappa_type) / 10)
  summary_data$type <- factor(summary_data$type, levels = paste("par", 1:15, sep = ""))

  pdf(paste(dir_plot, "pn_rho_pos_sigma", c("ind", "dep")[isigma], ".pdf", sep = ""), height = 7, width = 14)
  p <- summary_data %>% ggplot(aes(x = type, y = mean, type = factor(kappa))) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
    geom_errorbar(aes(ymin = lower, ymax = upper, col = kappa), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
    labs(x = "Parameter", y = "95% Credible Interval") +
    gg_theme +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "black", size = 1
    ) +
    scale_x_discrete(labels = c(
      "par1" = expression(rho[3][","][1]),
      "par2" = expression(rho[3][","][1]),
      "par3" = expression(rho[4][","][1]),
      "par4" = expression(rho[5][","][1]),
      "par5" = expression(rho[6][","][1]),
      "par6" = expression(rho[3][","][2]),
      "par7" = expression(rho[4][","][2]),
      "par8" = expression(rho[5][","][2]),
      "par9" = expression(rho[6][","][2]),
      "par10" = expression(rho[4][","][3]),
      "par11" = expression(rho[5][","][3]),
      "par12" = expression(rho[6][","][3]),
      "par13" = expression(rho[5][","][4]),
      "par14" = expression(rho[6][","][4]),
      "par15" = expression(rho[6][","][5])
    )) +
    labs(color = expression(kappa[j])) +
    facet_wrap(~n)
  print(p)
  dev.off()
}

# ========
# COPULA CASE
# ========



# load(paste(dir_data, "wc_simulations_results -  select_seed= 2 .Rdata", sep = ""))
# load(paste(dir_data, "pn_simulations_parameters -  select_seed= 1 .Rdata", sep = ""))
load(paste("simulations/output/cwc_simulations_results -  select_seed= 2 .Rdata", sep = ""))

n_vec <- c(20, 50, 100)
d_vec <- c(6)
w_n <- c()
w_d <- c()
w_k <- c()
for (isim in 1:length(out))
{
  if (out[[isim]]$d %in% d_vec) w_d <- c(w_d, isim)
  if (out[[isim]]$n %in% n_vec) w_n <- c(w_n, isim)
  if (out[[isim]]$rho[1] == out[[isim]]$rho[2]) w_k <- c(w_k, isim)
}
w_data <- intersect(w_n, w_d) %>% intersect(w_k)


## create dataset
app <- matrix(seq(1:(6^2)), ncol = 6)
index_cor <- app[lower.tri(app)]


h <- 1
isim <- w_data[1]

out_w <- out[[isim]]

kappa_type <- ifelse(out_w$rho[3] == 0.3, "1", ifelse(out_w$rho[3] == 0.6, "2", ifelse(out_w$rho[3] == out_w$rho[2], "3", "4")))

nsim_mcmc <- nrow(out_w$rho_out)
data_plot <- data.frame(
  d = out_w$d,
  n = out_w$n,
  # rho
  kappa_type = kappa_type,
  rho_sim = matrix(out_w$rho, ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
  rho = out_w$rho_out,
  # sigma
  cor_type = ifelse(out_w$Sigma_s[1, 2] == 0, "1", "2"),
  cor_sim = matrix(out_w$Sigma_s[lower.tri(app)], ncol = length(index_cor), nrow = nsim_mcmc, byrow = T),
  cor = out_w$sigma_s_out[, lower.tri(app)],
  # mu
  mu_type = "1",
  mu_sim = matrix(c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6)[1:out_w$d], ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
  mu = out_w$mu_out,
  index = h
)
for (isim in w_data)
{
  h <- h + 1
  out_w <- out[[isim]]

  kappa_type <- ifelse(out_w$rho[3] == 0.3, "1", ifelse(out_w$rho[3] == 0.6, "2", ifelse(out_w$rho[3] == out_w$rho[2], "3", "4")))

  nsim_mcmc <- nrow(out_w$rho_out)
  app_frame <- data.frame(
    d = out_w$d,
    n = out_w$n,
    # kappa
    kappa_type = kappa_type,
    rho_sim = matrix(out_w$rho, ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
    rho = out_w$rho_out,
    # sigma
    cor_type = ifelse(out_w$Sigma_s[1, 2] == 0, "1", "2"),
    cor_sim = matrix(out_w$Sigma_s[lower.tri(app)], ncol = length(index_cor), nrow = nsim_mcmc, byrow = T),
    cor = out_w$sigma_s_out[, lower.tri(app)],
    # mu
    mu_type = "1",
    mu_sim = matrix(c(0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6, 0, pi / 6, 2 * pi / 6, 3 * pi / 6, 4 * pi / 6, 5 * pi / 6)[1:out_w$d], ncol = out_w$d, nrow = nsim_mcmc, byrow = T),
    mu = out_w$mu_out,
    index = h
  )
  data_plot <- rbind(data_plot, app_frame)
}


# ========
# PLOTS
# ========


for (isigma in 1:2)
{
  #  isigma <- 1
  ## mu



  w_post <- which(colnames(data_plot) %in% paste("mu.", 1:out_w$d, sep = ""))
  w_sim <- which(colnames(data_plot) %in% paste("mu_sim.", 1:out_w$d, sep = ""))

  diff_var <- data_plot[, w_post] - data_plot[, w_sim]
  colnames(diff_var) <- paste("par", 1:out_w$d, sep = "")


  data_plot_app_1 <- data_plot %>%
    mutate(
      var = diff_var
    ) %>%
    unnest(
      cols = var,
      data = .
    ) %>%
    filter(.$cor_type == paste(isigma))

  data_plot_app <- data_plot_app_1 %>%
    select("n", "kappa_type", "cor_type", c(paste("par", 1:out_w$d, sep = ""))) %>%
    pivot_longer(
      cols = paste("par", 1:out_w$d, sep = ""),
      names_to = "type",
      values_to = "samp"
    )

  summary_data <- data_plot_app %>%
    group_by(type, cor_type, n, kappa_type) %>%
    summarize(
      mean = mean(ifelse(samp > pi, samp - 2 * pi, samp)),
      lower = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.025),
      upper = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.975)
    ) %>%
    mutate(lambda = case_when(
      kappa_type == "1" ~ "0.3",
      kappa_type == "2" ~ "0.6",
      kappa_type == "3" ~ "0.9"
    ))

  # summary_data <- summary_data %>%
  #   %>%
  #  mutate(kappa_order = as.numeric(kappa) + as.numeric(kappa_type) / 10)


  pdf(paste(dir_plot, "wc_mu_pos_sigma", c("ind", "dep")[isigma], ".pdf", sep = ""), height = 7, width = 14)
  p <- summary_data %>% ggplot(aes(x = type, y = mean, type = factor(lambda))) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
    geom_errorbar(aes(ymin = lower, ymax = upper, col = lambda), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
    labs(x = "Parameter", y = "95% Credible Interval") +
    gg_theme +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "black", size = 1
    ) +
    scale_x_discrete(labels = c(
      "par1" = expression(mu[1]),
      "par2" = expression(mu[2]),
      "par3" = expression(mu[3]),
      "par4" = expression(mu[4]),
      "par5" = expression(mu[5]),
      "par6" = expression(mu[6])
    )) +
    labs(color = expression(lambda[j])) +
    facet_wrap(~n)
  print(p)
  dev.off()



  ## kappa
  w_post <- which(colnames(data_plot) %in% paste("rho.", 1:out_w$d, sep = ""))
  w_sim <- which(colnames(data_plot) %in% paste("rho_sim.", 1:out_w$d, sep = ""))

  diff_var <- data_plot[, w_post] - data_plot[, w_sim]
  colnames(diff_var) <- paste("par", 1:out_w$d, sep = "")
  data_plot_app_1 <- data_plot %>%
    mutate(
      var = diff_var
    ) %>%
    unnest(
      cols = var,
      data = .
    ) %>%
    filter(.$cor_type == paste(isigma))

  data_plot_app <- data_plot_app_1 %>%
    select("n", "kappa_type", "cor_type", c(paste("par", 1:out_w$d, sep = ""))) %>%
    pivot_longer(
      cols = paste("par", 1:out_w$d, sep = ""),
      names_to = "type",
      values_to = "samp"
    )


  # data_plot_app <- cbind(data_plot_app_1, data_plot_app_2) %>% mutate(type = str_replace(type, "mu.", ""))


  summary_data <- data_plot_app %>%
    group_by(type, cor_type, n, kappa_type) %>%
    summarize(
      mean = mean(ifelse(samp > pi, samp - 2 * pi, samp)),
      lower = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.025),
      upper = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.975)
    ) %>%
    mutate(lambda = case_when(
      kappa_type == "1" ~ "0.3",
      kappa_type == "2" ~ "0.6",
      kappa_type == "3" ~ "0.9"
    ))


  pdf(paste(dir_plot, "wc_lambda_pos_sigma", c("ind", "dep")[isigma], ".pdf", sep = ""), height = 7, width = 14)
  p <- summary_data %>% ggplot(aes(x = type, y = mean, type = factor(lambda))) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
    geom_errorbar(aes(ymin = lower, ymax = upper, col = lambda), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
    labs(x = "Parameter", y = "95% Credible Interval") +
    gg_theme +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "black", size = 1
    ) +
    scale_x_discrete(labels = c(
      "par1" = expression(lambda[1]),
      "par2" = expression(lambda[2]),
      "par3" = expression(lambda[3]),
      "par4" = expression(lambda[4]),
      "par5" = expression(lambda[5]),
      "par6" = expression(lambda[6])
    )) +
    labs(color = expression(lambda[j])) +
    facet_wrap(~n)
  print(p)
  dev.off()

  ## rho
  w_post <- which(colnames(data_plot) %in% paste("cor.", 1:length(index_cor), sep = ""))
  w_sim <- which(colnames(data_plot) %in% paste("cor_sim.", 1:length(index_cor), sep = ""))

  diff_var <- data_plot[, w_post] - data_plot[, w_sim]
  colnames(diff_var) <- paste("par", 1:length(index_cor), sep = "")
  data_plot_app_1 <- data_plot %>%
    mutate(
      var = diff_var
    ) %>%
    unnest(
      cols = var,
      data = .
    ) %>%
    filter(.$cor_type == paste(isigma))

  data_plot_app <- data_plot_app_1 %>%
    select("n", "kappa_type", "cor_type", c(paste("par", 1:length(index_cor), sep = ""))) %>%
    pivot_longer(
      cols = paste("par", 1:length(index_cor), sep = ""),
      names_to = "type",
      values_to = "samp"
    )


  # data_plot_app <- cbind(data_plot_app_1, data_plot_app_2) %>% mutate(type = str_replace(type, "mu.", ""))


  summary_data <- data_plot_app %>%
    group_by(type, cor_type, n, kappa_type) %>%
    summarize(
      mean = mean(ifelse(samp > pi, samp - 2 * pi, samp)),
      lower = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.025),
      upper = quantile(ifelse(samp > pi, samp - 2 * pi, samp), 0.975)
    ) %>%
    mutate(lambda = case_when(
      kappa_type == "1" ~ "0.3",
      kappa_type == "2" ~ "0.6",
      kappa_type == "3" ~ "0.9"
    ))

  summary_data$type <- factor(summary_data$type, levels = paste("par", 1:15, sep = ""))


  # summary_data <- summary_data %>%
  #  mutate(kappa = case_when(
  #    kappa_type == "1" ~ "1",
  #    kappa_type == "2" ~ "10",
  #    kappa_type == "3" ~ "20",
  #    kappa_type == "4" & type == "par1" ~ "1",
  #    kappa_type == "4" & type == "par2" ~ "10",
  #    kappa_type == "4" & type == "par3" ~ "20",
  #    kappa_type == "4" & type == "par4" ~ "1",
  #    kappa_type == "4" & type == "par5" ~ "10",
  #    kappa_type == "4" & type == "par6" ~ "20"
  #  )) %>%
  #  mutate(kappa_order = as.numeric(kappa) + as.numeric(kappa_type) / 10)
  # summary_data %>% ggplot(aes(x = type, y = mean, type = factor(lambda))) +
  #  geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
  #  geom_errorbar(aes(ymin = lower, ymax = upper, col = lambda), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
  #  labs(x = "Parameter", y = "95% Credible Interval") +
  #  gg_theme +
  #  geom_hline(
  #    yintercept = 0, linetype = "dashed",
  #    color = "black", size = 1
  #  ) +
  #  scale_x_discrete(labels = c(
  #    "par1" = expression(rho[3][","][1]),
  #    "par2" = expression(rho[3][","][1]),
  #    "par3" = expression(rho[4][","][1]),
  #    "par4" = expression(rho[5][","][1]),
  #    "par5" = expression(rho[6][","][1]),
  #    "par6" = expression(rho[3][","][2]),
  #    "par7" = expression(rho[4][","][2]),
  #    "par8" = expression(rho[5][","][2]),
  #    "par9" = expression(rho[6][","][2]),
  #    "par10" = expression(rho[4][","][3]),
  #    "par11" = expression(rho[5][","][3]),
  #    "par12" = expression(rho[6][","][3]),
  #    "par13" = expression(rho[5][","][4]),
  #    "par14" = expression(rho[6][","][4]),
  #    "par15" = expression(rho[6][","][5])
  #  )) +
  #  facet_wrap(~n)

  pdf(paste(dir_plot, "wc_rho_pos_sigma", c("ind", "dep")[isigma], ".pdf", sep = ""), height = 7, width = 14)
  p <- summary_data %>% ggplot(aes(x = type, y = mean, type = factor(lambda))) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
    geom_errorbar(aes(ymin = lower, ymax = upper, col = lambda), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
    labs(x = "Parameter", y = "95% Credible Interval") +
    gg_theme +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "black", size = 1
    ) +
    scale_x_discrete(labels = c(
      "par1" = expression(rho[3][","][1]),
      "par2" = expression(rho[3][","][1]),
      "par3" = expression(rho[4][","][1]),
      "par4" = expression(rho[5][","][1]),
      "par5" = expression(rho[6][","][1]),
      "par6" = expression(rho[3][","][2]),
      "par7" = expression(rho[4][","][2]),
      "par8" = expression(rho[5][","][2]),
      "par9" = expression(rho[6][","][2]),
      "par10" = expression(rho[4][","][3]),
      "par11" = expression(rho[5][","][3]),
      "par12" = expression(rho[6][","][3]),
      "par13" = expression(rho[5][","][4]),
      "par14" = expression(rho[6][","][4]),
      "par15" = expression(rho[6][","][5])
    )) +
    labs(color = expression(lambda[j])) +
    facet_wrap(~n)
  print(p)
  dev.off()
}


# summary_data %>% ggplot(aes(x = type, y = mean, type = factor(lambda))) +
#  geom_point(size = 1.5, position = position_dodge(width = 0.5)) + # Plot mean
#  geom_errorbar(aes(ymin = lower, ymax = upper, col = lambda), width = 0.2, position = position_dodge(width = 0.5)) + # Plot credible interval
#  labs(x = "Parameter", y = "95% Credible Interval") +
#  gg_theme +
#  geom_hline(
#    yintercept = 0, linetype = "dashed",
#    color = "black", size = 1
#  ) +
#  scale_x_discrete(labels = c(
#    "par1" = expression(rho[3][","][1]),
#    "par2" = expression(rho[3][","][1]),
#    "par3" = expression(rho[4][","][1]),
#    "par4" = expression(rho[5][","][1]),
#    "par5" = expression(rho[6][","][1]),
#    "par6" = expression(rho[3][","][2]),
#    "par7" = expression(rho[4][","][2]),
#    "par8" = expression(rho[5][","][2]),
#    "par9" = expression(rho[6][","][2]),
#    "par10" = expression(rho[4][","][3]),
#    "par11" = expression(rho[5][","][3]),
#    "par12" = expression(rho[6][","][3]),
#    "par13" = expression(rho[5][","][4]),
#    "par14" = expression(rho[6][","][4]),
#    "par15" = expression(rho[6][","][5])
#  )) +
#  labs(color = expression(lambda[j]))
# facet_wrap(~n)

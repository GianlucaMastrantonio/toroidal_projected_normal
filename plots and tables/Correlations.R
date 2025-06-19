# circular_correlation <- function(theta, phi) {
#  # Ensure inputs are numeric vectors in radians
#  if (length(theta) != length(phi)) stop("Vectors must be of the same length")

#  n <- length(theta)

#  # Compute circular means
#  mean_theta <- atan2(mean(sin(theta)), mean(cos(theta)))
#  mean_phi <- atan2(mean(sin(phi)), mean(cos(phi)))

#  mean_theta <- 0
#  mean_phi <- 0

#  # Numerator
#  num <- sum(sin(theta - mean_theta) * sin(phi - mean_phi))

#  # Denominator
#  denom <- sqrt(
#    sum(sin(theta - mean_theta)^2) * sum(sin(phi - mean_phi)^2)
#  )

#  # Correlation
#  rho <- num / denom
#  return(rho)
# }
#################################################
# PLOTS OF CORRELATION COEFFICIENTS
# OF THE PROJECTED NORMAL DISTRIBUTION
#################################################
dir_plot <- "plots and tables/out/"
# library(cubature)
library(doParallel)
library(foreach)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

source("functions/general_functions.R")
gg_theme <- theme(
  plot.title = element_text(size = 18), # Increase title size
  axis.title = element_text(size = 22), # Increase axis title size
  axis.text = element_text(size = 20), # Increase axis text size
  legend.title = element_text(size = 25), # Increase legend title size
  legend.text = element_text(size = 16), # Increase legend text size
  legend.position = "bottom",
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  strip.text = element_text(size = 20)
)
########################
### TOriodal PN
######### ###############
# Monte Carlo
n_sim <- 10000
# Define a sequence of rho values
rho_seq <- seq(0, 0.99, length.out = 100) # originally 0.985 length.out=100 (changed for faster calculation)
rho_seq <- c(rho_seq, -rho_seq)
# Define kappa1 and kappa2 pairs to evaluate
kappa_pairs <- list(
  c(0, 0),
  c(0.49, 0.49),
  c(2.45, 2.45),
  c(0, 0.49),
  c(0, 2.45),
  c(0.49, 2.45)
)

# standard normales
x_gen_c <- mvrnorm(n = n_sim, mu = c(0, 0), Sigma = diag(1, 2))
x_gen_s <- mvrnorm(n = n_sim, mu = c(0, 0), Sigma = diag(1, 2))
corr_coef <- matrix(0, nrow = length(rho_seq), ncol = length(kappa_pairs))
for (i_rho in 1:length(rho_seq))
{
  for (i_kappa in 1:length(kappa_pairs))
  {
    rho <- rho_seq[i_rho]
    kappa <- kappa_pairs[[i_kappa]]
    # now i trasform the standard normals into normal i need
    Sigma_squared <- t(chol(matrix(c(1, rho, rho, 1), nrow = 2)))
    Sigma_squared_c <- abs(Sigma_squared)
    x_s <- t(Sigma_squared %*% t(x_gen_s))
    x_c <- t(Sigma_squared_c %*% t(x_gen_c))
    x_c[, 1] <- x_c[, 1] + kappa[1]
    x_c[, 2] <- x_c[, 2] + kappa[2]



    theta_1 <- atan2(x_s[, 1], x_c[, 1])
    theta_2 <- atan2(x_s[, 2], x_c[, 2])

    lambda <- (abs(mean(exp(1i * (theta_1 - theta_2)))) - abs(mean(exp(1i * (theta_1 + theta_2))))) / 2


    corr_coef[i_rho, i_kappa] <- lambda / max(c(mean(sin(theta_1)^2), mean(sin(theta_2)^2)))
    # corr_coef[i_rho, i_kappa] <- circular_correlation(theta_1, theta_2)
  }
}

data_plot <- data.frame(
  Correlation = c(corr_coef),
  rho = rep(rho_seq, times = length(kappa_pairs)),
  kappa = (rep(1:length(kappa_pairs), each = length(rho_seq)))
)

###

p1 <- data_plot %>%
  filter(kappa <= 3) %>%
  mutate(kappa = factor(kappa,
    levels = c(1, 2, 3),
    labels = c(
      "kappa[j]==''~kappa[k]==0",
      "kappa[j]==~''~kappa[k]==0.49",
      "kappa[j]==~''~kappa[k]==2.45"
      #  "kappa[1]==2.45~','~kappa[2]==2.45"
    )
  )) %>%
  ggplot(aes(x = rho, y = Correlation, color = kappa, linetype = kappa)) +
  geom_line(size = 2) +
  scale_color_discrete(labels = scales::parse_format()) +
  scale_linetype_discrete(labels = scales::parse_format()) +
  gg_theme +
  labs(x = expression(rho), y = "Circular Correlation", color = NULL, linetype = NULL) +
  guides(
    color = guide_legend(
      keywidth = 7, keyheight = 1.5, ,
      override.aes = list(linetype = c("solid", "dashed", "dotted"))
    )
  ) +
  theme(
    legend.direction = "vertical",
    legend.text = element_text(size = 20)
  ) +
  ylim(-1, 1)


pdf(paste(dir_plot, "Correlation1.pdf", sep = ""), height = 7, width = 7)
print(p1)
dev.off()



p2 <- data_plot %>%
  filter(kappa > 3) %>%
  mutate(kappa = factor(kappa,
    levels = c(4, 5, 6),
    labels = c(
      "kappa[j]==0~','~kappa[k]==0.49",
      "kappa[j]==0~','~kappa[k]==2.45",
      "kappa[j]==0.49~','~kappa[k]==2.45"
      #  "kappa[1]==2.45~','~kappa[2]==2.45"
    )
  )) %>%
  ggplot(aes(x = rho, y = Correlation, color = kappa, linetype = kappa)) +
  geom_line(size = 2) +
  scale_color_discrete(labels = scales::parse_format()) +
  scale_linetype_discrete(labels = scales::parse_format()) +
  gg_theme +
  labs(x = expression(rho), y = "Circular Correlation", color = NULL, linetype = NULL) +
  guides(
    color = guide_legend(
      keywidth = 7, keyheight = 1.5, ,
      override.aes = list(linetype = c("solid", "dashed", "dotted"))
    )
  ) +
  theme(
    legend.direction = "vertical",
    legend.text = element_text(size = 20)
  ) +
  ylim(-1, 1)

pdf(paste(dir_plot, "Correlation2.pdf", sep = ""), height = 7, width = 7)
print(p2)
dev.off()


########################
### Wrapped Cauchy
######### ###############
# Monte Carlo
n_sim <- 10000
# Define a sequence of rho values
rho_seq <- seq(0, 0.9999, length.out = 100) # originally 0.985 length.out=100 (changed for faster calculation)
rho_seq <- c(rho_seq, -rho_seq)
# Define kappa1 and kappa2 pairs to evaluate
wcpar_pairs <- list(
  c(0.0, 0.0),
  c(0.3, 0.3),
  c(0.9, 0.9),
  c(0.0, 0.3),
  c(0.0, 0.9),
  c(0.3, 0.9)
)

# standard normales
x_gen_c <- mvrnorm(n = n_sim, mu = c(0, 0), Sigma = diag(1, 2))
x_gen_s <- mvrnorm(n = n_sim, mu = c(0, 0), Sigma = diag(1, 2))
corr_coef <- matrix(0, nrow = length(rho_seq), ncol = length(kappa_pairs))
for (i_rho in 1:length(rho_seq))
{
  # print(i_rho)
  for (i_kappa in 1:length(kappa_pairs))
  {
    rho <- rho_seq[i_rho]
    wc_par <- wcpar_pairs[[i_kappa]]
    # now i trasform the standard normals into normal i need
    Sigma_squared <- t(chol(matrix(c(1, rho, rho, 1), nrow = 2)))
    Sigma_squared_c <- abs(Sigma_squared)
    x_s <- t(Sigma_squared %*% t(x_gen_s))
    x_c <- t(Sigma_squared_c %*% t(x_gen_c))




    theta_cop_1 <- atan2(x_s[, 1], x_c[, 1])
    theta_cop_2 <- atan2(x_s[, 2], x_c[, 2])
    theta_1 <- theta_cop_1
    theta_2 <- theta_cop_2
    for (iobs in 1:length(theta_cop_1))
    {
      theta_1[iobs] <- q_wc((theta_cop_1[iobs] %% (2 * pi)) / (2 * pi), 0, wc_par[1])
      theta_2[iobs] <- q_wc((theta_cop_2[iobs] %% (2 * pi)) / (2 * pi), 0, wc_par[2])
    }

    lambda <- (abs(mean(exp(1i * (theta_1 - theta_2)))) - abs(mean(exp(1i * (theta_1 + theta_2))))) / 2


    corr_coef[i_rho, i_kappa] <- circular_correlation(theta_1, theta_2)
  }
}

data_plot <- data.frame(
  Correlation = c(corr_coef),
  rho = rep(rho_seq, times = length(kappa_pairs)),
  kappa = (rep(1:length(kappa_pairs), each = length(rho_seq)))
)

###

p1 <- data_plot %>%
  filter(kappa <= 3) %>%
  mutate(kappa = factor(kappa,
    levels = c(1, 2, 3),
    labels = c(
      "lambda[j]==0~','~lambda[k]==0",
      "lambda[j]==0.3~','~lambda[k]==0.3",
      "lambda[j]==0.9~','~lambda[k]==0.9"
      #  "kappa[1]==2.45~','~kappa[2]==2.45"
    )
  )) %>%
  ggplot(aes(x = rho, y = Correlation, color = kappa, linetype = kappa)) +
  geom_line(size = 2) +
  scale_color_discrete(labels = scales::parse_format()) +
  scale_linetype_discrete(labels = scales::parse_format()) +
  gg_theme +
  labs(x = expression(rho), y = "Circular Correlation", color = NULL, linetype = NULL) +
  guides(
    color = guide_legend(
      keywidth = 7, keyheight = 1.5, ,
      override.aes = list(linetype = c("solid", "dashed", "dotted"))
    )
  ) +
  theme(
    legend.direction = "vertical",
    legend.text = element_text(size = 20)
  ) +
  ylim(-1, 1)


pdf(paste(dir_plot, "Correlation3.pdf", sep = ""), height = 7, width = 7)
print(p1)
dev.off()



p2 <- data_plot %>%
  filter(kappa > 3) %>%
  mutate(kappa = factor(kappa,
    levels = c(4, 5, 6),
    labels = c(
      "lambda[j]==0~','~lambda[k]==0.3",
      "lambda[j]==0~','~lambda[k]==0.9",
      "lambda[j]==0.3~','~lambda[k]==0.9"
      #  "kappa[1]==2.45~','~kappa[2]==2.45"
    )
  )) %>%
  ggplot(aes(x = rho, y = Correlation, color = kappa, linetype = kappa)) +
  geom_line(size = 2) +
  scale_color_discrete(labels = scales::parse_format()) +
  scale_linetype_discrete(labels = scales::parse_format()) +
  gg_theme +
  labs(x = expression(rho), y = "Circular Correlation", color = NULL, linetype = NULL) +
  guides(
    color = guide_legend(
      keywidth = 7, keyheight = 1.5, ,
      override.aes = list(linetype = c("solid", "dashed", "dotted"))
    )
  ) +
  theme(
    legend.direction = "vertical",
    legend.text = element_text(size = 20)
  ) +
  ylim(-1, 1)

pdf(paste(dir_plot, "Correlation4.pdf", sep = ""), height = 7, width = 7)
print(p2)
dev.off()

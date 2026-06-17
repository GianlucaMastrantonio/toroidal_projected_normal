library(ggplot2)
library(dplyr)
load("/Users/gianlucamastrantonio/Politecnico di Torino Staff Dropbox/Gianluca Mastrantonio/lavori/gitrepo/toroidal_projected_normal/simulations/output/post_analisi_results_cw.Rdata")


d_sel <- 1
n_sel <- 3
kappa_sel <- 1
select_sigma <- 1
w_data <- which(data_cw$select_d == d_sel & data_cw$select_n == n_sel & data_cw$select_kappa == kappa_sel & select_sigma == data_cw$select_sigma)
data_plot <- data_cw[w_data, ]
index_plot <- c()
sigma_mean <- c()
sigma_q1 <- c()
sigma_q3 <- c()
for (ii in 1:length(w_data))
{
  sigma_mean <- c(sigma_mean, list_res_cw[[w_data[ii]]]$sigma_mean)
  sigma_q1 <- c(sigma_q1, list_res_cw[[w_data[ii]]]$sigma_q1)
  sigma_q3 <- c(sigma_q3, list_res_cw[[w_data[ii]]]$sigma_q3)

  index_plot <- c(index_plot, rep(ii, length(list_res_cw[[w_data[ii]]]$sigma_mean)))
}

data_plot <- data.frame(index = index_plot, sigma_mean = unlist(sigma_mean), sigma_q1 = unlist(sigma_q1), sigma_q3 = unlist(sigma_q3))
library(ggplot2)
library(dplyr)
library(tidyr)

sigma_mat <- do.call(
  cbind,
  lapply(w_data, function(w) list_res_cw[[w]]$sigma_q3)
)

colnames(sigma_mat) <- paste0("index_", seq_along(w_data))

data_wide <- as.data.frame(sigma_mat)

pairs_plot <- combn(names(data_wide), 2, simplify = FALSE)

data_pairs <- bind_rows(
  lapply(pairs_plot, function(p) {
    data.frame(
      x = data_wide[[p[1]]],
      y = data_wide[[p[2]]],
      pair = paste(p[1], "vs", p[2])
    )
  })
)

ggplot(data_pairs, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~pair, scales = "free") +
  labs(
    x = "sigma_mean",
    y = "sigma_mean"
  ) +
  theme_bw()

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(mcclust.ext)



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
  #  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  strip.text = element_text(size = 20)
)


# data_plot <- data.frame(crps <- c(crps_tot[1, ], crps_tot[2, ]), k <- c(2:10, 2:10), model <- rep(c("PN", "WC"), each = 9))


# p1 <- data_plot %>% ggplot(aes(x = k, y = crps, type = factor(model), col = factor(model))) +
#  geom_line(size = 1.5) +
#  labs(x = "Parameter", y = "95% Credible Interval") +
#  gg_theme

# p1

### SECTION Maos

library(ggplot2)
library(maps)
library(ggplot2)
library(maps)
library(leaflet)
library(mapview)
library(webshot2)
if (!requireNamespace("osmdata", quietly = TRUE)) install.packages("osmdata")
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")

library(osmdata)
library(sf)
library(ggplot2)
library(maps)
if (!requireNamespace("leaflet", quietly = TRUE)) install.packages("leaflet")

load(file = "real data/data/gauge.Rdata")
library(leaflet)


circ_mean <- function(x) {
  atan2(mean(sin(x)), mean(cos(x)))
}

circular_cor_js <- function(theta_matrix) {
  theta_matrix <- as.matrix(theta_matrix)

  mu <- apply(theta_matrix, 2, circ_mean)

  centered_sin <- sin(sweep(theta_matrix, 2, mu, FUN = "-"))

  out <- cor(centered_sin)
  diag(out) <- 1

  out
}

theta <- theta_all
gauge_metadata <- gauge_metadata_all
corr_mat <- circular_cor_js(theta)

gauge_metadata$point_order <- seq_len(nrow(gauge_metadata))

gauge_metadata <- gauge_metadata_all
gauge_metadata <- gauge_metadata_all

gauge_metadata[6, "dec_long_va"] <- gauge_metadata[6, "dec_long_va"] + 0.05
gauge_metadata$point_order <- seq_len(nrow(gauge_metadata))


library(data.table)
library(leaflet)

thr <- 0.4

edges <- as.data.table(
  which(abs(corr_mat) > thr & upper.tri(corr_mat), arr.ind = TRUE)
)

setnames(edges, c("from", "to"))

edges[, correlation := corr_mat[cbind(from, to)]]

edges[, `:=`(
  lon_from = gauge_metadata$dec_long_va[from],
  lat_from = gauge_metadata$dec_lat_va[from],
  lon_to   = gauge_metadata$dec_long_va[to],
  lat_to   = gauge_metadata$dec_lat_va[to]
)]

edges[, line_width := 0.2 + 5 * abs(correlation)]

molt <- 0.7

m <- leaflet(gauge_metadata, width = 1200 * molt, height = 900 * molt) |>
  addProviderTiles(providers$Esri.WorldTopoMap)

for (ii in seq_len(nrow(edges))) {
  m <- m |>
    addPolylines(
      lng = c(edges$lon_from[ii], edges$lon_to[ii]),
      lat = c(edges$lat_from[ii], edges$lat_to[ii]),
      color = ifelse(edges$correlation[ii] > 0, "blue", "orange"),
      weight = edges$line_width[ii],
      opacity = 0.65,
      dashArray = ifelse(edges$correlation[ii] > 0, "", "18,20"),
      label = paste0("corr = ", round(edges$correlation[ii], 2))
    )
}

m <- m |>
  addCircleMarkers(
    lng = ~dec_long_va,
    lat = ~dec_lat_va,
    radius = 4,
    color = "red",
    fillColor = "red",
    fillOpacity = 0.85,
    stroke = TRUE,
    label = ~ paste0(column_index, ": ", site_no),
    popup = ~ paste0(
      "<b>Column:</b> ", column_index, "<br>",
      "<b>Gauge:</b> ", site_no, "<br>",
      "<b>Name:</b> ", station_nm
    )
  ) |>
  addLabelOnlyMarkers(
    lng = ~dec_long_va,
    lat = ~dec_lat_va,
    label = ~ as.character(point_order),
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "top",
      textOnly = TRUE,
      offset = c(0, -5),
      style = list(
        "color" = "black",
        "font-weight" = "bold",
        "font-size" = "15px",
        "background" = "transparent",
        "background-color" = "transparent",
        "border" = "none",
        "box-shadow" = "none"
      )
    )
  )
mapview::mapshot(
  m,
  file = "plots and tables/out/gauges_map.pdf",
  selfcontained = FALSE
)

## SECTION  correlation psoterior

load("real data/output_diagnostic/Ind=FALSE_do_small=1_ntry=40_molt_iter=1_model=tpnsetting_id=3.RData")
chains <- lapply(chain, function(ch) ch[["sigma_s_out"]])
x <- do.call(abind::abind, c(chains, along = 3))
x <- aperm(x, c(1, 3, 2))
dimnames(x) <- list(
  NULL,
  paste0("chain", 1:length(chain)),
  paste0("par", 1:dim(x)[3])
)

sigma_mean <- apply(x, c(3), mean)

rho_T <- function(rho) {
  sapply(rho, function(r) {
    if (!is.finite(r)) {
      return(NA_real_)
    }
    if (abs(r) >= 1) stop("rho must be in (-1, 1)")

    if (r == 0) {
      return(0)
    }

    z <- r^2

    K <- integrate(
      function(t) 1 / sqrt(1 - z * sin(t)^2),
      lower = 0,
      upper = pi / 2
    )$value

    E <- integrate(
      function(t) sqrt(1 - z * sin(t)^2),
      lower = 0,
      upper = pi / 2
    )$value

    val <- (E - (1 - z) * K)^2 / z

    sign(r) * val
  })
}
sigma_mean[sigma_mean >= 1] <- 0.999999
rho_mat <- rho_T(sigma_mean)


d <- 20
corr_long <- data.frame(var1 = rep(1:d, each = d), var2 = rep(1:d, times = d), corr = c(corr_mat))
library(ggplot2)


# data_plot <- data.frame(var = colMeans(sigma_s_out[, , k]), x = rep((1:d), each = d), y = rep((1:d), times = d))
p1 <- ggplot(corr_long, aes(x = var1, y = var2, fill = corr)) +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))

corr_long <- data.frame(var1 = rep(1:d, each = d), var2 = rep(1:d, times = d), corr = c(sigma_mean))
library(ggplot2)


# data_plot <- data.frame(var = colMeans(sigma_s_out[, , k]), x = rep((1:d), each = d), y = rep((1:d), times = d))
p2 <- ggplot(corr_long, aes(x = var1, y = var2, fill = corr)) +
  # geom_tile(color = "white") +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))


corr_long <- data.frame(var1 = rep(1:d, each = d), var2 = rep(1:d, times = d), corr = c(rho_mat))
library(ggplot2)


# data_plot <- data.frame(var = colMeans(sigma_s_out[, , k]), x = rep((1:d), each = d), y = rep((1:d), times = d))
# p3 <- ggplot(corr_long, aes(x = var1, y = var2, fill = corr)) +
#  # geom_tile(color = "white") +
#  geom_tile() +
#  scale_y_reverse() +
#  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))

p3 <- ggplot(corr_long, aes(x = factor(var1), y = factor(var2), fill = corr)) +
  geom_tile(color = "black", linewidth = 0.25) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(-1, 1),
    name = "corr"
  ) +
  coord_fixed() +
  theme_void() +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 14,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.text.y = element_text(
      color = "black",
      size = 14
    ),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.5, "cm")
  ) +
  labs(x = NULL, y = NULL)
# p3
# p1
# p2
# p3

pdf("plots and tables/out/corr_mat.pdf", height = 4.5 * 2, width = 4.5 * 2)
print(p3)
dev.off()



data_plot <- data.frame(corr = c(corr_mat)[c(corr_mat) < 1], rho = c(rho_mat)[c(corr_mat) < 1])
p1 <- ggplot(data_plot, aes(x = corr, y = rho)) +
  geom_point(size = 3) +
  gg_theme +
  labs(x = "Sample Correlation", y = "Posterior Mean of the Correlation") +
  gg_theme +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1)
pdf("plots and tables/out/corr_mat_points.pdf", height = 4.5 * 2, width = 4.5 * 2)
print(p1)
dev.off()

plot(c(corr_mat)[c(corr_mat) < 1], c(rho_mat)[c(corr_mat) < 1])
# plot(c(sigma_mean)[c(corr_mat) < 1], c(rho_mat)[c(corr_mat) < 1])
# plot(c(corr_mat), c(rho_mat))
# ggplot(corr_long, aes(x = var1, y = var2, fill = corr)) +
#  # geom_tile(color = "white") +
#  geom_tile() +
#  scale_y_reverse() +
#  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1))


# SECTION : DATA
#### select n
library(circular)
library(ggplot2)
library(data.table)
circular_density_all <- function(theta_all, n_grid = 512, kappa = 20, bw = 25) {
  theta_all <- as.matrix(theta_all)

  if (is.null(colnames(theta_all))) {
    colnames(theta_all) <- paste0("V", seq_len(ncol(theta_all)))
  }

  out <- vector("list", ncol(theta_all))

  for (j in seq_len(ncol(theta_all))) {
    x <- theta_all[, j]
    x <- x[is.finite(x)] %% (2 * pi)

    x_circ <- circular(
      x,
      units = "radians",
      modulo = "2pi",
      zero = 0,
      rotation = "counter"
    )

    dens <- density.circular(
      x_circ,
      n = n_grid,
      kernel = "vonmises",
      bw = 25,
      kappa = kappa
    )

    out[[j]] <- data.table::data.table(
      variable = colnames(theta_all)[j],
      theta = as.numeric(dens$x),
      density = as.numeric(dens$y)
    )
  }

  data.table::rbindlist(out)
}
theta <- theta_all
library(ggplot2)
library(tidyr)
library(dplyr)

theta_all_2pi <- theta_all %% (2 * pi)
dens_long <- circular_density_all(theta_all_2pi, n_grid = 1200, bw = 25, kappa = 15)

# dens_long <- circular_density_all(theta_all, n_grid = 512, bw = 25)

ggplot(dens_long, aes(x = theta, y = density)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~variable, scales = "free_y") +
  scale_x_continuous(
    breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
    labels = c("0", expression(pi / 2), expression(pi), expression(3 * pi / 2), expression(2 * pi))
  ) +
  theme_minimal() +
  labs(
    x = expression(theta),
    y = "Circular density"
  )


## here
n_bins <- 36


theta_plot <- theta
colnames(theta_plot) <- 1:20

theta_plot <- as.matrix(theta_plot)

# if (is.null(colnames(theta))) {
#  colnames(theta) <- paste0("V", seq_len(ncol(theta)))
# }

theta_long <- as.data.table(theta_plot)

theta_long[, row_id := seq_len(.N)]

theta_long <- melt(
  theta_long,
  id.vars = "row_id",
  variable.name = "variable",
  value.name = "theta"
)

theta_long[, theta := theta %% (2 * pi)]


p1 <- ggplot(theta_long, aes(x = theta)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(0, 2 * pi, length.out = n_bins + 1),
    fill = "grey70",
    color = "black"
  ) +
  facet_wrap(~variable, nrow = 4, ncol = 5, dir = "h") +
  scale_x_continuous(
    limits = c(0, 2 * pi),
    breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
    labels = c("0", expression(pi / 2), expression(pi), expression(3 * pi / 2), expression(2 * pi))
  ) +
  theme_minimal() +
  labs(
    x = expression(theta),
    y = "Density"
  )

pdf("plots and tables/out/hist_data.pdf", height = 7 * 1.2, width = 7 * 1)
print(p1)
dev.off()



# ! density estimates

chains <- lapply(chain, function(ch) ch[["sigma_s_out"]])
x <- do.call(abind::abind, c(chains, along = 3))
sigma_mcmc <- aperm(x, c(1, 3, 2))
dimnames(sigma_mcmc) <- list(
  NULL,
  paste0("chain", 1:length(chain)),
  paste0("par", 1:dim(sigma_mcmc)[3])
)


chains <- lapply(chain, function(ch) ch[["circ_mean"]])
x <- do.call(abind::abind, c(chains, along = 3))
mean_mcmc <- aperm(x, c(1, 3, 2))
dimnames(mean_mcmc) <- list(
  NULL,
  paste0("chain", 1:length(chain)),
  paste0("par", 1:dim(mean_mcmc)[3])
)


chains <- lapply(chain, function(ch) ch[["prec"]])
x <- do.call(abind::abind, c(chains, along = 3))
prec_mcmc <- aperm(x, c(1, 3, 2))
dimnames(prec_mcmc) <- list(
  NULL,
  paste0("chain", 1:length(chain)),
  paste0("par", 1:dim(prec_mcmc)[3])
)


library(circular)

lll <- 100
theta_seq <- circular(seq(0, 2 * pi, length.out = lll), units = "radians")

dens_est <- matrix(0, nrow = lll, ncol = d)
for (isim in 1:dim(prec_mcmc)[1])
{
  for (ic in 1:5)
  {
    for (id in 1:d)
    {
      # ss <- matrix(sigma_mcmc[isim, ic, ], nrow = d, ncol = d)
      Rmat <- matrix(0, nrow = 2, ncol = 2)
      Rmat[1, 1] <- cos(mean_mcmc[isim, ic, id])
      Rmat[2, 2] <- -sin(mean_mcmc[isim, ic, id])
      Rmat[1, 2] <- sin(mean_mcmc[isim, ic, id])
      Rmat[2, 1] <- cos(mean_mcmc[isim, ic, id])

      mu <- matrix(c(prec_mcmc[isim, ic, id], 0), nrow = 1) %*% Rmat
      dens_est[, id] <- dens_est[, id] + dpnorm(theta_seq, mu = mu, sigma = diag(1, 2))
    }
  }
}
dens_est <- dens_est / (dim(prec_mcmc)[1] * 5)

dens_grid <- seq(0, 2 * pi, length.out = nrow(dens_est))

dens_df <- as.data.frame(dens_est)

colnames(dens_df) <- levels(factor(theta_long$variable))

dens_long <- data.table::as.data.table(dens_df)
dens_long[, theta := dens_grid]

dens_long <- data.table::melt(
  dens_long,
  id.vars = "theta",
  variable.name = "variable",
  value.name = "density_est"
)

p1 <- ggplot(theta_long, aes(x = theta)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(0, 2 * pi, length.out = n_bins + 1),
    fill = "grey70",
    color = "black"
  ) +
  geom_line(
    data = dens_long,
    aes(x = theta, y = density_est),
    color = "red",
    linewidth = 0.8
  ) +
  facet_wrap(~variable, nrow = 4, ncol = 5, dir = "h") +
  scale_x_continuous(
    limits = c(0, 2 * pi),
    breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
    labels = c("0", expression(pi / 2), expression(pi), expression(3 * pi / 2), expression(2 * pi))
  ) +
  theme_minimal() +
  labs(
    x = expression(theta),
    y = "Density"
  )


pdf("plots and tables/out/hist_data.pdf", height = 7, width = 7 * 1.2)
print(p1)
dev.off()

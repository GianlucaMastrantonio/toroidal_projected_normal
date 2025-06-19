library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(gridExtra)


dir_plot <- "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/Applicazioni/Overleaf/Projected Normal Distributions on Torus/Figures/"


# ========
# PN
# ========

# cwc
load("real data/output/cwc.Rdata")
mean(c(crps_val))
mean(c(waic))
load("real data/output/tpn.Rdata")
mean(c(crps_val))
mean(c(waic))


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


## ========
## DATASET
## ========




p_list <- list()
theta_data_frame <- as.data.frame(theta)
for (icol in 1:ncol(theta_data_frame))
{
  theta_data_frame[, icol] <- ifelse(theta_data_frame[, icol] > pi, theta_data_frame[, icol] - 2 * pi, theta_data_frame[, icol])
}
# Ensure the columns are named correctly


colnames(theta_data_frame) <- c("Staz_PTAT2_Wind", "Staz_42035_Wind", "Staz_42019_Wind", "Staz_42019_Wave", "Staz_42020_Wind", "Staz_42020_Wave")[c(1, 2, 3, 5, 4, 6)]
names_col <- colnames(theta_data_frame)
names_plot <- c("PTAT2 - Wind", "42035 - Wind", "42019 - Wind", "42019 - Wave", "42020 - Wind", "42020 - Wave")[c(1, 2, 3, 5, 4, 6)]
h <- 1
for (i in 1:5)
{
  for (j in (i + 1):6)
  {
    p_list[[h]] <- theta_data_frame %>% ggplot(aes_string(x = names_col[i], y = names_col[j])) +
      geom_point(size = 1) +
      gg_theme +
      labs(
        x = names_plot[i],
        y = names_plot[j]
      ) +
      scale_x_continuous(limits = c(-pi, pi)) + # Set the x-axis to be between 0 and 2*pi
      scale_y_continuous(limits = c(-pi, pi)) +
      theme(
        legend.position = "none", # Removes the legend
        legend.title = element_blank() # Optional: to remove the legend title
      ) +
      coord_fixed(ratio = 1)
    h <- h + 1
  }
}

molt_fig <- 2
pdf(paste(dir_plot, "data_all.pdf", sep = ""), height = 4.5 * 2 * molt_fig, width = 3 * 2 * molt_fig)
print(grid.arrange(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], p_list[[5]], p_list[[6]], p_list[[7]], p_list[[8]], p_list[[9]], p_list[[10]], p_list[[11]], p_list[[12]], p_list[[13]], p_list[[14]], p_list[[15]], ncol = 3, nrow = 5))
dev.off()
# Plot using ggplot2


##### Parameters
id <- 1
data_plot <- data.frame(var = c(mu_out[, id]), iter = 1:dim(mu_out)[1], type = 1)
for (id in 1:d)
{
  data_plot <- rbind(data_plot, data.frame(var = c(mu_out[, id]), iter = 1:dim(mu_out)[1], type = id))
}
p1 <- data_plot %>% ggplot(aes(x = var)) +
  geom_density(adjust = 2, key_glyph = draw_key_path, size = 1) +
  gg_theme +
  scale_linetype_manual(values = c("1" = "solid", "2" = "dashed", "3" = "dotted", "4" = "dotdash", "5" = "twodash")) +
  guides(linetype = guide_legend(keywidth = 4, keyheight = 2)) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_sqrt() +
  facet_wrap(~type, labeller = labeller(type = c("1" = names_plot[1], "2" = names_plot[2], "3" = names_plot[3], "4" = names_plot[4], "5" = names_plot[5], "6" = names_plot[6])))

pdf(paste(dir_plot, "post_par1.pdf", sep = ""), height = 7 * 0.9 * 1.8, width = 7 * 1 * 1.8)
print(p1)
dev.off()


id <- 1
data_plot <- data.frame(var = c(rho_out[, id]), iter = 1:dim(mu_out)[1], type = 1)
for (id in 1:d)
{
  data_plot <- rbind(data_plot, data.frame(var = c(rho_out[, id]), iter = 1:dim(mu_out)[1], type = id))
}
p1 <- data_plot %>% ggplot(aes(x = var)) +
  geom_density(adjust = 2, key_glyph = draw_key_path, size = 1) +
  gg_theme +
  scale_linetype_manual(values = c("1" = "solid", "2" = "dashed", "3" = "dotted", "4" = "dotdash", "5" = "twodash")) +
  guides(linetype = guide_legend(keywidth = 4, keyheight = 2)) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_sqrt() +
  facet_wrap(~type, labeller = labeller(type = c("1" = names_plot[1], "2" = names_plot[2], "3" = names_plot[3], "4" = names_plot[4], "5" = names_plot[5], "6" = names_plot[6])))

pdf(paste(dir_plot, "post_par2.pdf", sep = ""), height = 7 * 0.9 * 1.8, width = 7 * 1 * 1.8)
print(p1)
dev.off()





data_plot <- data.frame(var = colMeans(sigma_s_out[, ]), x = rep((1:d), each = d), y = rev(rep((1:d), times = d)), type = 1)


p1 <- data_plot %>% ggplot(aes(x = x, y = y, fill = var)) +
  geom_tile() +
  gg_theme +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5, limits = c(0, 1)) +
  theme(legend.position = "bottom") +
  guides(linetypeinc = guide_legend(keywidth = 4, keyheight = 2)) +
  theme(
    legend.key.height = unit(1, "cm"), # Adjust height of legend items
    legend.key.width = unit(3, "cm"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_x_continuous(
    breaks = 1:6, # Set the x-axis ticks at specific values
    labels = names_plot # Custom labels for the ticks
  ) +
  scale_y_continuous(
    breaks = 6:1, # Set the x-axis ticks at specific values
    labels = names_plot # Custom labels for the ticks
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1) # Rotate y-axis labels by 90 degrees
  )
p1

pdf(paste(dir_plot, "post_var2.pdf", sep = ""), height = 7 * 0.9 * 1.8, width = 7 * 1 * 1.8)
print(p1)
dev.off()



# data_plot <- data.frame(z = zeta_map, time = 1:nrow(theta))

# p1 <- data_plot %>% ggplot(aes(x = time, y = z)) +
#  geom_line(size = 2) +
#  gg_theme +
#  theme(legend.position = "bottom") +
#  guides(linetypeinc = guide_legend(keywidth = 4, keyheight = 2)) +
#  theme(
#    legend.key.height = unit(1, "cm"), # Adjust height of legend items
#    legend.key.width = unit(3, "cm"),
#    legend.title = element_blank(),
#    axis.title.x = element_blank(),
#    axis.title.y = element_blank()
#  )
# p1
# pdf(paste("zeta.pdf", sep = ""), height = 7 * 0.9, width = 7 * 1 * 2)
# print(p1)
# dev.off()

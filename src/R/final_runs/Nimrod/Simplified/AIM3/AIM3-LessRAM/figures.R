# Final figures for honours report
# 0) Conceptual diagram of the model space
# 1) A) Prove we are at equilibrium with variance over time vs selection strength
# 1) B) Define our space of models between HOC and Gaussian with end-sim Va vs tau
# 2) Diagram of what our models are (adapted from Walsh and Lynch figure 28.1), the space between the discrete types
# 3) One aspect of these models calculations are the effects as Ne -> Inf, we summarise that with delmu, how does delmu influence mean variance at equilibrium with locisigma
# 4) Why does variance matter? Consequences for adaptation around an optimum, distance vs delmu and locisigma
# 5) Unpack variance more: distributions of additive effects at equilibrium - deleterious mutation and locisigma
# 6) Conceptual diagram of distances from optimum and random scatter according to rates of delmu at low and high locisigma

# Set constants etc. import libraries


source("../../AIM1/AIM-1_R/src_G_mat.R")
source("../../AIM1/AIM-1_R/src_plot.R")

# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

library(plyr)
library(tidyverse)


# Import data

d_eucdist_c <- readRDS("d_eucdist_c.RDS")

d_raw_c <- readRDS("d_raw_c.RDS")

d_raw_end <- readRDS("d_raw_end.RDS")



# Axis labels

pr_lab <- "Rate of pleiotropy"
pc_lab <- "Mutational pleiotropic correlation"
r_lab <- "Recombination rate"
d_lab <- "Deleterious mutation rate"
ls_lab <- "Additive effect size"
t_lab <- "Selection strength (\u03C4)"
dist_lab <- "\u03B4\u0305"
var_lab <- "\u03C3\u0305\u00B2"
cov_lab <- "Cov(X, Y)"


# 1) A) Prove we are at equilibrium in respect to variance/covariance

# Mean data across seeds

d_mean_vartime <- d_raw_c[, c(1, 3, 5:6, 15:26)] %>%
  group_by(gen, modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error, var = var))

plot_vartime <- ggplot(d_mean_vartime, aes(x = gen, y = varmean_groupmean)) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  scale_color_npg() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = var_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# 1) B) Define the space of models with selection strength vs mutation rate balance
# Plot null vs selection separately


dplot_var_t <- d_raw_end[d_raw_end$tau != 0, c(15:16, 20)] %>%
  group_by(tau) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))

dplot_var_t[dplot_var_t$tau == 0,]$tau <- 10000


plot_var_t <- ggplot(dplot_var_t, aes(x = tau, y = varmean_mean)) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ x, col = "#03A9F4", fill = "#5FBBE6") +
  scale_y_continuous(limits = c(0, 800)) +
  theme_classic() +
  ggtitle(t_lab) +
  labs(x = t_lab, y = var_lab) +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_t


# 3) Mean variance at equilibrum vs delmu and locisigma

dplot_var_cont <- d_raw_end[, c(3, 5:6, 15:26)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(mean = mean, se = std.error, var = var))

plot_var_cont.d.ls <- ggplot(dplot_var_cont, aes(x = delmu, y = varmean_mean)) +
  facet_grid(locisigma.cat~.) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  #  ggtitle(d_lab) +
  labs(x = d_lab, y = var_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_cont.d.ls


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_var_cont.d.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_var_cont.d.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_var_cont.d.ls)


# 4) Mean distance at equilibrium vs delmu and locisigma

dplot_eucdist_cont <- d_eucdist_fingen[, c(3:16)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))

plot_dist_cont.d.ls <- ggplot(dplot_eucdist_cont, aes(x = delmu, y = dist_mean)) +
  facet_grid(locisigma.cat~.) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  #  ggtitle(d_lab) +
  labs(x = d_lab, y = dist_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_dist_cont.d.ls


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_dist_cont.d.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_dist_cont.d.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_dist_cont.d.ls)




# Final figures for honours report
# -1) Latin Square sampling
# 0) Conceptual diagram of the model space
# 1) A) Prove we are at equilibrium with variance over time vs selection strength - also distance over time
# 1) B) Define our space of models between HOC and Gaussian with end-sim Va vs tau
# 2) Diagram of what our models are (adapted from Walsh and Lynch figure 28.1), the space between the discrete 
# types: sampling different levels of mutation vs deleterious mutation
# 3) Additive variance is supposed to predict how well you are able to get to an optimum. How does the additive
# effect size distribution influence trait variance? Tau, rwide?
# 4) Why does variance matter? Consequences for adaptation around an optimum, distance vs locisigma
# 5) Unpack variance more: distributions of additive effects at equilibrium - locisigma, tau, rwide (?)
# 6) Conceptual diagram of distances from optimum and random scatter according to rates of delmu at low and high locisigma
# 7) Percentage of models that reach the optimum


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

d_mean_eucdist <- readRDS("d_mean_eucdist.RDS")

d_mean_var <- readRDS("d_mean_var.RDS")

d_raw_end <- readRDS("d_raw_end.RDS")

# Use only from delmu = 0.25 to 0.75 to ensure we've got polygenic architecture: delmu = 1 means 50% fewer mutations
# contributing to trait than delmu = 0

d_eucdist_c <- d_eucdist_c[d_eucdist_c$delmu <= 0.75 & d_eucdist_c$delmu >= 0.25,]

d_eucdist_fingen <- d_eucdist_fingen[d_eucdist_fingen$delmu <= 0.75 & d_eucdist_fingen$delmu >= 0.25,]

d_raw_c <- d_raw_c[d_raw_c$delmu <= 0.75 & d_raw_c$delmu >= 0.25,]

d_raw_end <- d_raw_end[d_raw_end$delmu <= 0.75 & d_raw_end$delmu >= 0.25,]

levels(d_eucdist_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_eucdist_fingen$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_end$tau.cat) <- c("Null", "Strong", "Medium", "Weak")


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
po_lab <- "P(o)"

# Colours

cs <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1, length.out = 100))
cs <- c("black", cs[c(100, 15, 1)])


# 1) A) Prove we are at equilibrium in respect to variance/covariance
# Mean across parameters

dplot_varmean_time <- d_raw_c[, c(1, 15:16, 26)] %>%
  group_by(gen, tau.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))

plot_vartime <- ggplot(dplot_varmean_time, aes(x = gen, y = varmean_mean, colour = tau.cat, fill = tau.cat)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = varmean_mean - 1.96*varmean_se,
    ymax = varmean_mean + 1.96*varmean_se
  ), linetype = 0, alpha = 0.2) +
  scale_colour_manual(values = cs) +
  ggtitle("A") + 
  scale_fill_manual(values = cs) +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = var_lab, colour = t_lab, fill = t_lab) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))


plot_vartime  

# 1) A.5) Euc dist over time

dplot_dist_time <- d_eucdist_c[, c(1, 4, 13)] %>%
  group_by(gen, tau.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))



plot_disttime_t <- ggplot(dplot_dist_time, aes(x = gen, y = dist_mean, col = tau.cat, fill = tau.cat)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = dist_mean - 1.96*dist_se,
    ymax = dist_mean + 1.96*dist_se
  ), linetype = 0, alpha = 0.2) +
  scale_colour_manual(values = cs) +
  scale_fill_manual(values = cs) +
  ggtitle("B") + 
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = dist_lab, col = t_lab, fill = t_lab) +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))


plot_disttime_t


# Plot Fig 1 and 1.5 together

library(patchwork)

plot_dist_vartime <-  plot_vartime / plot_disttime_t

ggsave(filename = "dist+vartime.t.png", plot = plot_dist_vartime, width = 8, height = 12, dpi = 800)


# 1) B) Define the space of models with selection strength vs mutation rate balance
# Plot null vs selection separately


dplot_var_t <- d_raw_end[, c(15:16, 20)] %>%
  group_by(tau) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))

dplot_var_t[dplot_var_t$tau == 0,]$tau <- 1500


plot_var_t <- ggplot(d_raw_end, aes(x = tau.cat, y = varmean)) +
  geom_point() +
  geom_boxplot()+
  scale_colour_manual(values = cs) +
  theme_classic() +
  labs(x = t_lab, y = var_lab) +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_t

ggsave(filename = "var.t.png", plot = plot_var_t, width = 8, height = 8, dpi = 800)


plot_cov_t <- ggplot(d_raw_end, aes(x = tau.cat, y = covmean)) +
  geom_point() +
  geom_boxplot()+
  geom_smooth(se = T, method = "lm", formula = y ~ x, col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  ggtitle(t_lab) +
  labs(x = t_lab, y = cov_lab) +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_cov_t




# 3) A) Mean variance at equilibrium vs locisigma B) euclidean distance vs locisigma

dplot_var_cont <- d_raw_end[, c(3, 5:6, 15:26)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(mean = mean, se = std.error, var = var))

# 3)A: Variance


# LS

plot_var_cont.ls <- ggplot(dplot_var_cont, aes(x = locisigma, y = varmean_mean)) +
  geom_point(data = d_raw_end, mapping = aes(x=locisigma, y = varmean), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  ggtitle("A") +
  labs(x = ls_lab, y = var_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))

plot_var_cont.ls



# r

plot_var_cont.r <- ggplot(dplot_var_cont, aes(x = rwide, y = varmean_mean)) +
  geom_point(data = d_raw_end, mapping = aes(x=rwide, y = varmean), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  ggtitle("B") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 0.00012, 0.00004), labels = seq(0, 1.2, 0.4)) +
  labs(x = expression(bold(Recombination~rate~(x*"10"^"-4"))), y = var_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))

plot_var_cont.r

# t

plot_var_cont.t <- ggplot(dplot_var_cont[dplot_var_cont$tau.cat != "Null",], aes(x = tau, y = varmean_mean)) +
  geom_point(data = d_raw_end[d_raw_end$tau.cat != "Null",], mapping = aes(x=tau, y = varmean), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  labs(x = t_lab, y = var_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))

plot_var_cont.t


var_ls.r <- plot_var_cont.ls + plot_var_cont.r

ggsave(filename = "var_ls.r.png", plot = var_ls.r, width = 12, height = 8, dpi = 800)

#############################################################
# 3)B: Euclidean Distance

dplot_eucdist_cont <- d_eucdist_fingen[, c(3:16)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))

# LS

plot_dist_cont.ls <- ggplot(dplot_eucdist_cont, aes(x = locisigma, y = dist_mean)) +
  geom_point(data = d_eucdist_fingen, mapping = aes(x=locisigma, y = distance), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  ggtitle("C") +
  labs(x = ls_lab, y = dist_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))

plot_dist_cont.ls


# r

plot_dist_cont.r <- ggplot(dplot_eucdist_cont, aes(x = rwide, y = dist_mean)) +
  geom_point(data = d_eucdist_fingen, mapping = aes(x=rwide, y = distance), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  scale_x_continuous(breaks = seq(0, 0.00012, 0.00004), labels = seq(0, 1.2, 0.4)) +
  theme_classic() +
  ggtitle("D") +
  labs(x = expression(bold(Recombination~rate~(x*"10"^"-4"))), y = dist_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))

plot_dist_cont.r



eucdist_ls.r <- plot_dist_cont.ls | plot_dist_cont.r

ggsave(filename = "eucdist_ls.r.png", plot = eucdist_ls.r, width = 12, height = 8, dpi = 800)

var_dist_ls.r <- ( plot_var_cont.ls | plot_var_cont.r ) / ( plot_dist_cont.ls | plot_dist_cont.r ) 

ggsave(filename = "var_dist_ls.r.png", plot = var_dist_ls.r, width = 12, height = 12, dpi = 800)


##############
# Covariance
# LS

plot_cov_cont.ls <- ggplot(dplot_var_cont, aes(x = locisigma, y = covmean_mean)) +
  geom_point(data = d_raw_end, mapping = aes(x=locisigma, y = covmean), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  #  ggtitle(d_lab) +
  labs(x = ls_lab, y = cov_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_cov_cont.ls



# r

plot_cov_cont.r <- ggplot(dplot_var_cont, aes(x = rwide, y = covmean_mean)) +
  geom_point(data = d_raw_end, mapping = aes(x=rwide, y = covmean), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  #  ggtitle(d_lab) +
  labs(x = r_lab, y = cov_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_cont.r


# 7) Percentage of models that reach the optimum

dplot_po <- d_eucdist_fingen[,c(4:5, 9:10, 13:15)] %>%
  group_by(delmu, locisigma, delmu.cat, tau, locisigma.cat, tau.cat) %>%
  summarise_all(list(po = percent_dist), tol=1) 

d_eucdist_fingen$atopt <- d_eucdist_fingen$distance
d_eucdist_fingen[d_eucdist_fingen$distance <= 1.0,]$atopt <- 1
d_eucdist_fingen[d_eucdist_fingen$distance > 1.0,]$atopt <- 0


# Combine medium and high locisigma (same trend)
levels(dplot_po$locisigma.cat) <- c("Small", "Large", "Large")

plot_po.t.ls <- ggplot(dplot_po, aes(x = locisigma, y = po)) +
  facet_grid(tau.cat~.) +
  geom_point(data = d_eucdist_fingen, mapping = aes(x=locisigma, y = atopt), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 3), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  #  ggtitle(d_lab) +
  labs(x = ls_lab, y = po_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 18, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.x = unit(2, "lines"))

plot_po.t.ls



# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add labels

library(gtable)
library(grid)

# Labels 
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_po.t.ls)

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
plot_gtab_po.t.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_po.t.ls)

ggsave(filename = "po.t.ls.png", plot = plot_gtab_po.t.ls, width = 6, height = 8, dpi = 800)


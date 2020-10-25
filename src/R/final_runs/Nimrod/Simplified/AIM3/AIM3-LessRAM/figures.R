# Final figures for honours report
# -1) Latin Square sampling
# 0) Conceptual diagram of the model space
# 1) A) Prove we are at equilibrium with variance over time vs selection strength - also distance over time
# 1) B) Define our space of models between HOC and Gaussian with end-sim Va vs tau
# 2) Diagram of what our models are (adapted from Walsh and Lynch figure 28.1), the space between the discrete types
# 3) One aspect of these models calculations are the effects as Ne -> Inf, we summarise that with delmu, how does delmu influence mean variance at equilibrium with locisigma
# 4) Why does variance matter? Consequences for adaptation around an optimum, distance vs delmu and locisigma
# 5) Unpack variance more: distributions of additive effects at equilibrium - deleterious mutation and locisigma
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
  scale_fill_manual(values = cs) +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = var_lab, colour = t_lab, fill = t_lab) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
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
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = dist_lab, col = t_lab, fill = t_lab) +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))


plot_disttime_t


# Plot Fig 1 and 1.5 together

library(patchwork)

plot_dist_vartime <- plot_disttime_t / plot_vartime

ggsave(filename = "dist+vartime.t.png", plot = plot_dist_vartime, width = 8, height = 12, dpi = 800)


# 1) B) Define the space of models with selection strength vs mutation rate balance
# Plot null vs selection separately


dplot_var_t <- d_raw_end[, c(15:16, 20)] %>%
  group_by(tau) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))

dplot_var_t[dplot_var_t$tau == 0,]$tau <- 1500


plot_var_t <- ggplot(d_raw_c[d_raw_c$gen == 150000,], aes(x = tau, y = varmean, group = 1)) +
  geom_point() +
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


plot_cov_t <- ggplot(dplot_var_t, aes(x = tau, y = covmean_mean, group = 1)) +
  geom_point() +
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




# 3) A) Mean variance at equilibrium vs delmu and locisigma under neutral drift B) Mean variance at equilibrum vs delmu and locisigma under selection

dplot_var_cont <- d_raw_end[, c(3, 5:6, 15:26)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(mean = mean, se = std.error, var = var))

# 3)A 

dplot_var_cont_null <- dplot_var_cont[dplot_var_cont$tau.cat == "Null",]

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


ggsave(filename = "var_d.ls.png", plot = plot_gtab_var_cont.d.ls, width = 12, height = 8, dpi = 800)

#############################################################
# 3)B: Selection model

dplot_var_cont_sel <- dplot_var_cont[dplot_var_cont$tau.cat != "Null",]


plot_var_cont.d.ls.s <- ggplot(dplot_var_cont_sel, aes(x = delmu, y = varmean_mean)) +
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

plot_var_cont.d.ls.s


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_var_cont.d.ls.s)

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
plot_gtab_var_cont.d.ls.s <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_var_cont.d.ls.s)



##############
# Covariance

plot_cov_cont.d.ls <- ggplot(dplot_var_cont, aes(x = delmu, y = covmean_mean)) +
  facet_grid(locisigma.cat~.) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  #  ggtitle(d_lab) +
  labs(x = d_lab, y = cov_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_cov_cont.d.ls


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_cov_cont.d.ls)

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
plot_gtab_cov_cont.d.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_cov_cont.d.ls)

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

ggsave(filename = "dist_d.ls.png", plot = plot_gtab_dist_cont.d.ls, width = 12, height = 8, dpi = 800)



# 7) Percentage of models that reach the optimum

dplot_po <- d_eucdist_fingen[d_eucdist_fingen$tau.cat != "Null" ,c(4:5, 9:10, 13:15)] %>%
  group_by(delmu, locisigma, delmu.cat, tau, locisigma.cat, tau.cat) %>%
  summarise_all(list(po = percent_dist), tol=1) 

# Combine medium and high locisigma (same trend)
levels(dplot_po$locisigma.cat) <- c("Low", "High", "High")

plot_po.d.ls <- ggplot(dplot_po, aes(x = delmu, y = po)) +
  facet_grid(locisigma.cat~tau.cat) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  #  ggtitle(d_lab) +
  labs(x = d_lab, y = po_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 18, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.x = unit(2, "lines"))

plot_po.d.ls



# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add labels

library(gtable)
library(grid)

# Labels 
labelT = t_lab
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_po.d.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab_po.d.ls <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab_po.d.ls)

ggsave(filename = "po.d.ls.png", plot = plot_gtab_po.d.ls, width = 12, height = 8, dpi = 800)


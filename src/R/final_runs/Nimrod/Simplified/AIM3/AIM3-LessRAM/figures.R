# Final figures for honours report
# Supplementary:
# S1) Latin Square sampling: (A) Null (B) Sel
# S2) Heterozygosity
# S3) Validation of mu as mutation rate rather than deleterious mutation rate
# Table S1) Parameter sampling ranges and description of parameters

# Main text:
# 1) Conceptual diagram of the model space (HOC/Gaussian etc.) (Intro)
# 2) A) Prove we are at equilibrium with variance over time vs HOC vs Gaussian
# 2) B) Same as above except distance over time
# 3) A) Additive variance is supposed to predict how well you are able to get to an optimum. How does the additive
# effect size distribution influence trait variance, and how does it differ between models?
# 3) B) Genetic covariance isn't well studied in this regard: do the models behave differently with trait covariance?
# 4) Why does variance matter? Consequences for adaptation around an optimum, distance vs locisigma with model
# 5) Unpack variance more: distributions of additive effects at equilibrium - locisigma, model: this also tells us if the mutational variance is higher than standing genetic variance
# 6) Conceptual diagram of distances from optimum and random scatter according to different models at low and high locisigma
# 7) What lies between? Plot "Other" models between HOC and Gaussian with end-sim Va vs model type/locisigma

# Tables
# 1) ANOVA tables for Va, cov(x,y), sqrt(distance)
# 2) ANOVA tables for kurtosis, mean, variance of allele frequencies

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

levels(d_eucdist_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_eucdist_fingen$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_raw_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_raw_end$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")



levels(d_eucdist_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_eucdist_fingen$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_end$tau.cat) <- c("Null", "Strong", "Medium", "Weak")

# Comparing between models: Add a new column which is the interaction between selection strength and delmu (actually just mu)

d_eucdist_c$COA.cat <- interaction(d_eucdist_c$delmu.cat, d_eucdist_c$tau.cat)

levels(d_eucdist_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                 "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                 "Other", "Other")

d_eucdist_fingen <- d_eucdist_c[d_eucdist_c$gen == 150000,]

                    
d_raw_end$COA.cat <- interaction(d_raw_end$delmu.cat, d_raw_end$tau.cat)

levels(d_raw_end$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                               "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                               "Other", "Other")

d_raw_c$COA.cat <- interaction(d_raw_c$delmu.cat, d_raw_c$tau.cat)

levels(d_raw_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                             "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                             "Other", "Other")

d_mean_eucdist$COA.cat <- interaction(d_mean_eucdist$delmu.cat, d_mean_eucdist$tau.cat)

levels(d_mean_eucdist$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                               "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                               "Other", "Other")

d_mean_var$COA.cat <- interaction(d_mean_var$delmu.cat, d_mean_var$tau.cat)

levels(d_mean_var$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                               "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                               "Other", "Other")


# Axis labels

pr_lab <- "Rate of pleiotropy"
pc_lab <- "Mutational pleiotropic correlation"
r_lab <- "Recombination rate"
d_lab <- "Mutation rate"
ls_lab <- "Additive effect size distribution (\u03B1)"
t_lab <- "Selection strength (\u03C4)"
dist_lab <- "\u03B4\u0305"
var_lab <- "\u03C3\u0305\u00B2"
cov_lab <- "Cov(X, Y)"
po_lab <- "P(o)"
m_lab <- "Allelic effect model"

# Colours

cs <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1, length.out = 100))
cs <- c("black", cs[c(100, 15, 1)])


# 1) A) Prove we are at equilibrium in respect to variance/covariance
# Mean across parameters

dplot_varmean_time <- d_raw_c[, c(1, 15:16, 27)] %>%
  group_by(gen, COA.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))

plot_vartime <- ggplot(dplot_varmean_time[dplot_varmean_time$COA.cat != "Other",], aes(x = gen, y = varmean_mean, colour = COA.cat, fill = COA.cat)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = varmean_mean - 1.96*varmean_se,
    ymax = varmean_mean + 1.96*varmean_se
  ), linetype = 0, alpha = 0.2) +
  scale_colour_manual(values = c("black", "blue", "red"))+
  ggtitle("A") + 
  scale_fill_manual(values = c("black", "blue", "red")) +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = var_lab, colour = m_lab, fill = m_lab) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))


plot_vartime  

# 1) A.5) Euc dist over time

dplot_dist_time <- d_eucdist_c[, c(1, 4, 17)] %>%
  group_by(gen, COA.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))



plot_disttime_t <- ggplot(dplot_dist_time[dplot_dist_time$COA.cat != "Other",], aes(x = gen, y = dist_mean, col = COA.cat, fill = COA.cat)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = dist_mean - 1.96*dist_se,
    ymax = dist_mean + 1.96*dist_se
  ), linetype = 0, alpha = 0.2) +
  scale_colour_manual(values = c("black", "blue", "red")) +
  scale_fill_manual(values = c("black", "blue", "red")) +
  ggtitle("B") + 
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = dist_lab, col = m_lab, fill = m_lab) +
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


dplot_var_t <- d_raw_end[, c(15:16, 27)] %>%
  group_by(COA.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))

dplot_var_t[dplot_var_t$tau == 0,]$tau <- 1500


plot_var_t <- ggplot(d_raw_end[d_raw_end$COA.cat != "Other",], aes(x = COA.cat, y = varmean)) +
  geom_point() +
  geom_boxplot()+
  scale_colour_manual(values = c("black", "blue", "red")) +
  theme_classic() +
  ggtitle("A") +
  labs(x = m_lab, y = var_lab) +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))

plot_var_t


plot_cov_t <- ggplot(d_raw_end[d_raw_end$COA.cat != "Other",], aes(x = COA.cat, y = covmean)) +
  geom_point() +
  geom_boxplot()+
  geom_smooth(se = T, method = "lm", formula = y ~ x, col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  ggtitle("B") +
  labs(x = m_lab, y = cov_lab) +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))

plot_cov_t


plot_dist_varcov <-  plot_var_t + plot_cov_t

ggsave(filename = "dist+varend.t.png", plot = plot_dist_varcov, width = 12, height = 8, dpi = 800)


# 3) A) Mean variance at equilibrium vs locisigm and model B) cov vs locisigma and model

dplot_var_cont <- d_raw_end[, c(3, 5:6, 15:27)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat, COA.cat) %>%
  summarise_all(list(mean = mean, se = std.error, var = var))

# 3)A: Variance


# LS

plot_var_cont.ls <- ggplot(dplot_var_cont[dplot_var_cont$COA.cat != "Other",], aes(x = locisigma, y = varmean_mean)) +
  facet_grid(COA.cat~.) +
  geom_point(data = d_raw_end[d_raw_end$COA.cat != "Other",], mapping = aes(x=locisigma, y = varmean), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  labs(x = ls_lab, y = var_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_var_cont.ls


library(gtable)
library(grid)

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelR = "Allelic effect model"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_var_cont.ls)

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
plot_gtab_var_cont.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))


# Draw it
grid.newpage()
grid.draw(plot_gtab_var_cont.ls)


# Covariance

plot_cov_cont.ls <- ggplot(dplot_var_cont[dplot_var_cont$COA.cat != "Other",], aes(x = locisigma, y = covmean_mean)) +
  facet_grid(COA.cat~.) +
  geom_point(data = d_raw_end[d_raw_end$COA.cat != "Other",], mapping = aes(x=locisigma, y = covmean), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  labs(x = ls_lab, y = cov_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 8),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_cov_cont.ls



library(gtable)
library(grid)

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelR = "Allelic effect model"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_cov_cont.ls)

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
plot_gtab_cov_cont.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))



# Add corner labels before joining: thanks to: https://stackoverflow.com/a/29863172/13586824
plot_gtab_var_cont.ls <- arrangeGrob(plot_gtab_var_cont.ls, top = textGrob("A", 
                                                  x = unit(0, "npc"), 
                                                  y = unit(1, "npc"), 
                                                  just=c("left","top"), 
                                                  gp = gpar(col = "black", fontsize=22, fontface = "bold")))

plot_gtab_cov_cont.ls <- arrangeGrob(plot_gtab_cov_cont.ls, top = textGrob("B", 
                                                                           x = unit(0, "npc"), 
                                                                           y = unit(1, "npc"), 
                                                                           just=c("left","top"), 
                                                                           gp = gpar(col = "black", fontsize=22, fontface = "bold")))


varcov_ls.m <- grid.arrange(plot_gtab_var_cont.ls, plot_gtab_cov_cont.ls, nrow = 2)

ggsave(filename = "varcov_ls.m.png", plot = varcov_ls.m, width = 8, height = 16, dpi = 800)

#############################################################
# 4): Euclidean Distance

dplot_eucdist_cont <- d_eucdist_fingen[, c(3:17)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat, COA.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error, dist_var = var))

# LS

plot_dist.ls.m <- ggplot(dplot_eucdist_cont[dplot_eucdist_cont$COA.cat != "Other",], aes(x = locisigma, y = dist_mean)) +
  facet_grid(COA.cat~.) +
  geom_point(data = d_eucdist_fingen[d_eucdist_fingen$COA.cat != "Other",], mapping = aes(x=locisigma, y = distance), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2), col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  labs(x = ls_lab, y = dist_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

  plot_dist.ls.m

  
  # Thanks to: https://stackoverflow.com/a/37292665/13586824
  # Add rwide label
  # Labels 
  labelR = "Allelic effect model"
  
  # Get the ggplot grob
  plot_gtab <- ggplotGrob(plot_dist.ls.m)
  
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
  plot_gtab_dist.ls.m <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
  
  


ggsave(filename = "plot_dist.ls.m.png", plot = plot_gtab_dist.ls.m, width = 8, height = 12, dpi = 800)


##############

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



# Supplementary: S3 Latin Hypercube Sampling
library(GGally) # for multivariate scatterplots
lscombos_null <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM1/lscombos_null.csv")
lscombos_null <- lscombos_null[,-1]
lscombos_sel <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")
lscombos_sel <- lscombos_sel[,-1]

lscombos_null$rwide <- scales::rescale(lscombos_null$rwide)
lscombos_null$locisigma <- scales::rescale(lscombos_null$locisigma)
lscombos_null$pleiorate <- scales::rescale(lscombos_null$pleiorate)
lscombos_null$delmu <- scales::rescale(lscombos_null$delmu)
lscombos_null$pleiocov <- scales::rescale(lscombos_null$pleiocov)



names(lscombos_null) <- c(rwide = "Recombination\nrate (per locus)", 
                          locisigma = "Additive effect\nsize (\u03B1)", 
                          pleiorate = "Rate of\npleiotropy", 
                          delmu = "Deleterious\nmutation rate", 
                          pleiocov = "Mutational\ncorrelation")



plot_lhc_null <- ggpairs(lscombos_null, xlab = "Scaled parameter value", ylab = "Scaled parameter value",
                         lower = list(continuous = wrap("points", size = 0.1))) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  labs(tag = "A") +
  theme(text = element_text(size = 22, face = "bold"),
        panel.spacing = unit(1, "lines")
)
plot_lhc_null

ggsave("plot_lhc_null.png", plot_lhc_null, height = 12, width = 12, dpi = 400)



# Selection

lscombos_sel$rwide <- scales::rescale(lscombos_sel$rwide)
lscombos_sel$locisigma <- scales::rescale(lscombos_sel$locisigma)
lscombos_sel$pleiorate <- scales::rescale(lscombos_sel$pleiorate)
lscombos_sel$delmu <- scales::rescale(lscombos_sel$delmu)
lscombos_sel$pleiocov <- scales::rescale(lscombos_sel$pleiocov)
lscombos_sel$tau <- scales::rescale(lscombos_sel$tau)



names(lscombos_sel) <- c(rwide = "Recombination\nrate (per locus)", 
                          locisigma = "Additive effect\nsize (\u03B1)", 
                          pleiorate = "Rate of\npleiotropy", 
                          delmu = "Deleterious\nmutation rate", 
                          pleiocov = "Mutational\ncorrelation",
                          tau = "Selection\nstrength (\u03C4)")



plot_lhc_sel <- ggpairs(lscombos_sel, xlab = "Scaled parameter value", ylab = "Scaled parameter value",
                         lower = list(continuous = wrap("points", size = 0.1))) +
  theme_classic() +
  labs(tag = "B") +
  scale_x_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  theme(text = element_text(size = 22, face = "bold"),
        panel.spacing = unit(1, "lines")
  )
plot_lhc_sel

ggsave("plot_lhc_sel.png", plot_lhc_sel, height = 16, width = 16, dpi = 300)

library(patchwork)
plot_lhc <- grid.arrange(plot_lhc_null, plot_lhc_sel)

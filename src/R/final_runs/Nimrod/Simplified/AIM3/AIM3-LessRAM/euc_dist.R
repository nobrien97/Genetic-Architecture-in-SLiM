# Euclidean distances from optimum at final generation
#################################################################


# Distance from optimum - Euclidean distance in 8D space

source("../../AIM1/AIM-1_R/src_G_mat.R")


# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

# Figure labels
pr_lab <- "Rate of pleiotropy"
pc_lab <- "Mutational pleiotropic correlation"
r_lab <- "Recombination rate"
d_lab <- "Deleterious mutation rate"
ls_lab <- "Additive effect size"
t_lab <- "Selection strength (\u03C4)"

# Import data - optima and means etc.

d_sel <- read.csv("d_sel.csv")

d_opt <- data.table::fread("F:/Uni/AIM3/OUTPUT/out_8T_stabsel_opt_c.csv", header = F, integer64="character")

# Linux version

d_opt <- data.table::fread("/mnt/f/Uni/AIM3/OUTPUT/out_8T_stabsel_opt_c.csv", header = F, integer64="character")


names(d_opt) <- c("seed", "modelindex", paste0("opt", 0:7))
d_opt$seed <- as.numeric(d_opt$seed) # For sorting

ls_popmeans <- MCmean_gen(d_sel, d_sel$modelindex, 4)
ls_opt <- opt_gen(d_opt, d_opt$modelindex)

ls_eucdists <- MCeuc_dist(d_sel, d_opt, 4)



# Convert to data frame for plotting: need to take outer list as gen, next level as seed, and next level as model

library(tidyverse)

d_eucdist <- data.frame(
  gen = rep(unique(d_sel$gen), each = length(unique(d_sel$seed))*length(unique(d_sel$modelindex))),
  seed = sort(rep(unique(d_sel$seed), each = length(unique(d_sel$modelindex)))),
  modelindex = sort(unique(d_sel$modelindex)),
  distance = unlist(ls_eucdists)
)

# Verify it has sorted correctly
eg_means <- d_sel[,39:46][d_sel$seed == 2767825865 & d_sel$modelindex == 60,]
eg_opts <- d_opt[,3:10][d_opt$seed == 2767825865 & d_opt$modelindex == 60,]

names(eg_means) <- names(eg_opts)

dist(rbind(eg_means, eg_opts))
#> dist(rbind(eg_means, eg_opts)) eg_means <- d_sel[,39:46][d_sel$seed == 16601004 & d_sel$modelindex == 5,] eg_opts <- d_opt[,3:10][d_opt$seed == 16601004 & d_opt$modelindex == 5,]
#454
#1 208.2571

# Correct, we'll do another one just to double check

#> dist(rbind(eg_means, eg_opts)) eg_means <- d_sel[,39:46][d_sel$seed == 2767825865 & d_sel$modelindex == 60,] eg_opts <- d_opt[,3:10][d_opt$seed == 2767825865 & d_opt$modelindex == 60,]

#5936
#1 152.2622

# Also correct, so the data frame is ready to be filled with predictors
# Load in ls_cpombos if we haven't already
# Linux version
ls_combos_sel <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")


d_eucdist$delmu <- rep(ls_combos_sel$delmu, times = 100)
d_eucdist$rwide <- rep(ls_combos_sel$rwide, times = 100)
d_eucdist$pleiocov <- rep(ls_combos_sel$pleiocov, times = 100) # Repeat 100 times, for each seed
d_eucdist$pleiorate <- rep(ls_combos_sel$pleiorate, times = 100)
d_eucdist$locisigma <- rep(ls_combos_sel$locisigma, times = 100)
d_eucdist$tau <- rep(ls_combos_sel$tau, times = 100)

# Split into low/medium/high groups

d_eucdist$delmu.cat <- cut(d_eucdist$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$rwide.cat <- cut(d_eucdist$rwide, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$pleiocov.cat <- cut(d_eucdist$pleiocov, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$pleiorate.cat <- cut(d_eucdist$pleiorate, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$locisigma.cat <- cut(d_eucdist$locisigma, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$tau.cat <- cut(d_eucdist$tau, breaks = 3, labels = c("High", "Medium", "Low"))

# Reorder so predictors are first

d_eucdist <- d_eucdist[,c(1:3, 5:16, 4)]

write.csv(d_eucdist, "d_eucdist.csv", row.names = F)

d_eucdist <- read.csv("d_eucdist.csv")

# Order factors

d_eucdist$delmu.cat <- factor(d_eucdist$delmu.cat, levels = c("Low", "Medium", "High"))
d_eucdist$rwide.cat <- factor(d_eucdist$rwide.cat, levels = c("Low", "Medium", "High"))
d_eucdist$pleiocov.cat <- factor(d_eucdist$pleiocov.cat, levels = c("Low", "Medium", "High"))
d_eucdist$pleiorate.cat <- factor(d_eucdist$pleiorate.cat, levels = c("Low", "Medium", "High"))
d_eucdist$locisigma.cat <- factor(d_eucdist$locisigma.cat, levels = c("Low", "Medium", "High"))
d_eucdist$tau.cat <- factor(d_eucdist$tau.cat, levels = c("High", "Medium", "Low"))


# Type III lm/ANOVA followed by least squares means contrasts as post hoc


library(emmeans)
library(estimatr)
lm_eucdist <- lm_robust(distance ~ delmu.cat*rwide.cat*locisigma.cat + tau.cat,
                        data = d_eucdist)

summary(lm_eucdist)

emm_dist_contr_d.r.t <- pairs(pairs(emmeans(lm_eucdist, ~ delmu.cat * rwide.cat | tau.cat,
                                      at = list(delmu.cat = c("Low", "High"),
                                                rwide.cat = c("Low", "High"),
                                                tau.cat = c("Low")))), by = NULL)

emm_dist_d.r.ls <- emmeans(lm_eucdist, pairwise ~ delmu.cat * rwide.cat * locisigma.cat)
emm_dist_d.t <- emmeans(lm_eucdist, pairwise ~ delmu.cat | tau.cat)


pairs(pairs(emmeans(m, ~ f2|f1)), by = NULL)

emm_dist_pr.pc.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * pleiocov.cat | tau.cat)
emm_dist_pr.r.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * rwide.cat | tau.cat)

emm_dist_contr_pr.r.t <- pairs(pairs(emmeans(lm_eucdist, ~ pleiorate.cat * rwide.cat | tau.cat,
                                            at = list(pleiorate.cat = c("Low", "High"),
                                                      rwide.cat = c("Low", "High"),
                                                      tau.cat = c("High")))), by = NULL)


emm_dist_pr.d.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * delmu.cat | tau.cat)


emm_dist_contr_pr.d.t <- pairs(pairs(emmeans(lm_eucdist, ~ pleiorate.cat * delmu.cat | tau.cat,
                                             at = list(pleiorate.cat = c("Low", "High"),
                                                       delmu.cat = c("Low", "High"),
                                                       tau.cat = c("Low", "High")))), by = NULL)



emm_dist_pc.d.t <- emmeans(lm_eucdist, pairwise ~ pleiocov.cat * delmu.cat | tau.cat)
emm_dist_pr.ls.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * locisigma.cat | tau.cat)

emm_dist_contr_pr.ls.t <- pairs(pairs(emmeans(lm_eucdist, ~ pleiorate.cat * locisigma.cat | tau.cat,
                                             at = list(pleiorate.cat = c("Low", "High"),
                                                       locisigma.cat = c("Low", "High"),
                                                       tau.cat = c("High")))), by = NULL)



emm_dist_pc.ls.t <- emmeans(lm_eucdist, pairwise ~ pleiocov.cat * locisigma.cat | tau.cat)
emm_dist_pc.r.t <- emmeans(lm_eucdist, pairwise ~ pleiocov.cat * rwide.cat | tau.cat)


############################################################################################################
############################################################################################################
############################################################################################################
# Delmu * rwide * ls
dplot_eucdist_d.r.ls <- d_eucdist[,c(10, 11, 14, 16)] %>%
  group_by(delmu.cat, rwide.cat, locisigma.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_r.d.ls <- ggplot(dplot_eucdist_d.r.ls, aes(x = delmu.cat, y = dist_mean, fill = rwide.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Euclidean distance from optimum", fill = r_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_r.d.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.r.d.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.r.d.ls)


############################################################################################################
############################################################################################################
############################################################################################################

# Delmu * rwide * ls

plot_eucdist_d.r.ls <- ggplot(dplot_eucdist_d.r.ls, aes(x = rwide.cat, y = dist_mean, fill = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = r_lab, y = "Euclidean distance from optimum", fill = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_d.r.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.d.r.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.d.r.ls)


############################################################################################################
############################################################################################################
############################################################################################################





# Plot euclidean distance from optimum: separate figure for each parameter, coloured lines for tau bin

########################################
# Colour scale for graph points
cs <- scales::seq_gradient_pal("blue", "#ff8c00", "Lab")(seq(0,1, length.out = 100))
cs <- cs[c(1, 15, 45, 100)]

# Or use npg
library(ggsci)


# Delmu * rwide
dplot_eucdist_delmu.rwide <- d_eucdist[,c(10, 11, 15:16)] %>%
  group_by(delmu.cat, rwide.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_d.r <- ggplot(dplot_eucdist_delmu.rwide, aes(x = delmu.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(rwide.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = r_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_d.r)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.d.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.d.r)


# # # # # # # # # # # # # # # # # 


# pleiorate * pleiocov
dplot_eucdist_pleiorate.pleiocov <- d_eucdist[,c(12, 13, 15:16)] %>%
  group_by(pleiorate.cat, pleiocov.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.pc <- ggplot(dplot_eucdist_pleiorate.pleiocov, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(pleiocov.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = pc_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.pc)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.pc <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.pc)


# # # # # # # # # # # # # # # # # 


# pleiorate * rwide
dplot_eucdist_pleiorate.rwide <- d_eucdist[,c(11, 13, 15:16)] %>%
  group_by(pleiorate.cat, rwide.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.r <- ggplot(dplot_eucdist_pleiorate.rwide, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(rwide.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = r_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.r)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.r)


# # # # # # # # # # # # # # # # # 

# pleiocov * rwide
dplot_eucdist_pleiocov.rwide <- d_eucdist[,c(11, 12, 15:16)] %>%
  group_by(pleiocov.cat, rwide.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pc.r <- ggplot(dplot_eucdist_pleiocov.rwide, aes(x = pleiocov.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(rwide.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pc_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiocov label
# Labels 
library(grid)
library(gtable)
labelR = r_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pc.r)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pc.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pc.r)

# # # # # # # # # # # # # # # # # 


# pleiorate * delmu
dplot_eucdist_pleiorate.delmu <- d_eucdist[,c(10, 13, 15:16)] %>%
  group_by(pleiorate.cat, delmu.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.d <- ggplot(dplot_eucdist_pleiorate.delmu, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(delmu.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = d_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.d)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.d)


# # # # # # # # # # # # # # # # # 


# pleiocov * delmu
dplot_eucdist_pleiocov.delmu <- d_eucdist[,c(10, 12, 15:16)] %>%
  group_by(pleiocov.cat, delmu.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pc.d <- ggplot(dplot_eucdist_pleiocov.delmu, aes(x = pleiocov.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(delmu.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pc_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiocov label
# Labels 
library(grid)
library(gtable)
labelR = d_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pc.d)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pc.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pc.d)



# # # # # # # # # # # # # # # # # 


# pleiorate * locisigma
dplot_eucdist_pleiorate.locisigma <- d_eucdist[,c(13, 14, 15:16)] %>%
  group_by(pleiorate.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.ls <- ggplot(dplot_eucdist_pleiorate.locisigma, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.ls)


# # # # # # # # # # # # # # # # # 


# pleiocov * locisigma
dplot_eucdist_pleiocov.locisigma <- d_eucdist[,c(12, 14, 15:16)] %>%
  group_by(pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pc.ls <- ggplot(dplot_eucdist_pleiocov.locisigma, aes(x = pleiocov.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pc_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiocov label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pc.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pc.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pc.ls)

########################################################################
########################################################################
########################################################################

# Plot the whole lot together
library(gridExtra)

grid.arrange(plot_gtab_eucdist.pr.r, plot_gtab_eucdist.pc.r, newpage = T)
grid.arrange(plot_gtab_eucdist.pr.d, plot_gtab_eucdist.pc.d, newpage = T)
grid.arrange(plot_gtab_eucdist.pr.ls, plot_gtab_eucdist.pc.ls, newpage = T)
grid.arrange(plot_gtab_eucdist.d.r, plot_gtab_eucdist.pr.pc, newpage = T)



########################################################################################
# Over time: need to get pop means, mean var and mean cov so file shouldn't be too big
# Will need to do with selection model

# For null, will have to calculate the optimum after the fact
# phenomeans at generation 50,000
# optimum = phenomeans+(phenomeans*(3.1746) # mu*nloci*50,000 = 6.3e-6*100*50,000 = 31.5
# 

# Test here first with a small section, do big null file on supercomputer

d_null_256 <- data.table::fread("/mnt/f/Uni/AIM1/OUTPUT/out_8T_null_means_256.csv", header = F, integer64="character")


names(d_null_256)[1:6] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu")
# d_null_256$seed <- as.factor(d_null_256$seed)

names(d_null_256)[7:34] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_null_256)[35:42] <- paste0("mean", 0:7)

names(d_null_256)[43:50] <- paste0("var", 0:7)

names(d_null_256)[51:78] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_null_256)[79:106] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_null_256)[107] <- "H"

d_null_opt <- d_null_256[d_null_256$gen == 50000]
d_null_opt <- cbind(d_null_opt[, c(2:3)], d_null_opt[, c(35:42)] + (d_null_opt[, c(35:42)]*(100/31.5)))
names(d_null_opt) <- c("seed", "modelindex", paste0("opt", 0:7))


d_null$varmean <- rowMeans(d_null[, c(43:50)])
d_null$covmean <- rowMeans(d_null[, c(51:78)])

d_null <- d_null[, -c(43:107)]
write.csv(d_null_opt, "d_null_opt.csv")
write.csv(d_null, "d_null_time.csv")


# Linux version

d_sel <- data.table::fread("/mnt/f/Uni/AIM3/OUTPUT/out_8T_stabsel_means_c2.csv", header = F, integer64="character", fill = T, select = c(1:107))


# Remove duplicate rows
library(tidyverse)
d_sel <- d_sel %>% distinct()


names(d_sel)[1:7] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu", "tau")
# d_null$seed <- as.factor(d_null$seed)

names(d_sel)[8:35] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_sel)[36:43] <- paste0("mean", 0:7)

names(d_sel)[44:51] <- paste0("var", 0:7)

names(d_sel)[52:79] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_sel)[80:107] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

d_sel$seed <- as.numeric(d_sel$seed)

d_sel$varmean <- rowMeans(d_sel[, c(44:51)])
d_sel$covmean <- rowMeans(d_sel[, c(52:79)])

d_sel <- d_sel[, -c(8:35, 44:107)]

d_sel <- d_sel[order(d_sel$modelindex),]

# Add actual pleiocov line (value from the latin hypercube)
#ls_combos <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/pilot_runs/Pilot_Project/lscombos_sel.csv")

# Linux version
ls_combos_sel <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")

d_sel$locisigma <- rep(ls_combos_sel$locisigma, each = 20100)


write.csv(d_sel, "d_sel_time.csv", row.names = F)

d_sel <- read.csv("d_sel_time.csv")
d_sel <- d_sel[, -1]

d_sel_test <- arrange(d_sel, gen, modelindex, seed)
d_opt_test <- arrange(d_opt, modelindex, seed)

d_sel_test <- as.numeric(d_sel_test[1, c(8:15)])
d_opt_test <- as.numeric(d_opt_test[1, c(3:10)])
dist(rbind(d_sel_test, d_opt_test))


# sel distances
d_sel_opt <- data.table::fread("/mnt/f/Uni/AIM3/OUTPUT/out_8T_stabsel_opt_c.csv", header = F, integer64="character")
names(d_sel_opt) <- c("seed", "modelindex", paste0("opt", 0:7))
d_sel_opt$seed <- as.numeric(d_sel_opt$seed)


##########################################################################################
# Run on supercomputer: extract means, and optima, do distance

d_eucdist_c <- readRDS("d_eucdist_c.RDS")
d_mean_eucdist <- readRDS("d_mean_eucdist.RDS")
d_mean_eucdist_notau_null <- readRDS("d_mean_eucdist_notau_null.RDS")
d_mean_eucdist_notau_sel <- readRDS("d_mean_eucdist_notau_sel.RDS")


# Linear model: euclidean distance vs parameters at final time point
d_eucdist_fingen <- d_eucdist_c[d_eucdist_c$gen == 150000,]

library(emmeans)
library(estimatr)

eucdist_fingen_lm <- lm_robust(distance ~ delmu.cat * rwide.cat * locisigma.cat + 
                                 tau.cat * delmu.cat + tau.cat * rwide.cat + tau.cat * locisigma.cat,
data = d_eucdist_fingen)

summary(eucdist_fingen_lm)

emm_dist_contr_d.r.t <- pairs(pairs(emmeans(lm_eucdist, ~ delmu.cat * rwide.cat | tau.cat,
                                            at = list(delmu.cat = c("Low", "High"),
                                                      rwide.cat = c("Low", "High"),
                                                      tau.cat = c("Low")))), by = NULL)

emm_dist_d.r.ls <- emmeans(lm_eucdist, pairwise ~ delmu.cat * rwide.cat * locisigma.cat)
emm_dist_d.t <- emmeans(lm_eucdist, pairwise ~ delmu.cat | tau.cat)


# Plot euclidean distance from optimum: separate figure for each parameter, coloured lines for tau bin
# Bars at final timepoint, then a line graph of distance over time

########################################
# Colour scale for graph points
cs <- scales::seq_gradient_pal("blue", "#ff8c00", "Lab")(seq(0,1, length.out = 100))
cs <- cs[c(1, 15, 45, 100)]

d_eucdist_nullbar <- d_mean_eucdist_notau_null[d_mean_eucdist_notau_null$gen == 150000,]
d_eucdist_selbar <- d_mean_eucdist_notau_sel[d_mean_eucdist_notau_sel$gen == 150000,]

# Everything together
d_mean_eucdist_end <- d_mean_eucdist[d_mean_eucdist$gen == 150000,]


# Or use npg
library(ggsci)


# Null models and sel models: null models should have little difference in distance between treatments
# Maybe smaller effects mean they don't go as far, but delmu shouldn't impact it nor should recom

# Sel models decrease distance with increasing delmu, recombination reverses this effect somewhat
# Under low delmu, recombination can slightly reduce distance to optimum

plot_eucdist_r.d.ls.n <- ggplot(d_eucdist_nullbar, aes(x = delmu.cat, y = dist_mean, fill = rwide.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Euclidean distance from optimum", fill = r_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_r.d.ls.n)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.r.d.ls.n <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.r.d.ls.n)



plot_eucdist_r.d.ls.s <- ggplot(d_eucdist_selbar, aes(x = delmu.cat, y = dist_mean, fill = rwide.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Euclidean distance from optimum", fill = r_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_r.d.ls.s)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.r.d.ls.n <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.r.d.ls.s)



plot_eucdist_d.r.ls.n <- ggplot(d_eucdist_nullbar, aes(x = rwide.cat, y = dist_mean, fill = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = r_lab, y = "Euclidean distance from optimum", fill = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_d.r.ls.n)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.d.r.ls.n <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.d.r.ls.n)




plot_eucdist_d.r.ls.s <- ggplot(d_eucdist_selbar, aes(x = rwide.cat, y = dist_mean, fill = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = r_lab, y = "Euclidean distance from optimum", fill = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_d.r.ls.s)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.d.r.ls.s <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.d.r.ls.s)


# Distance over time

d_mean_eucdist_notau_null <- readRDS("d_mean_eucdist_notau_null.RDS")
d_mean_eucdist_notau_sel <- readRDS("d_mean_eucdist_notau_sel.RDS")

plot_disttime_d.ls.n <- ggplot(d_mean_eucdist_notau_null, aes(x = gen, y = dist_mean, col = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_line(position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Euclidean distance from optimum", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_disttime_d.ls.n)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_disttime.d.ls.n <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_disttime.d.ls.n)


# Average across all selection strengths
plot_disttime_d.ls.s <- ggplot(d_mean_eucdist_notau_sel, aes(x = gen, y = dist_mean, col = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_line(position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Euclidean distance from optimum", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_disttime_d.ls.s)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_disttime.d.ls.s <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_disttime.d.ls.s)

ggsave(filename = "disttime_d.ls.av_sel.png", plot = plot_gtab_disttime.d.ls.s, width = 12, height = 8, dpi = 800)





# Check if tau matters

d_mean_eucdist_lowmed <- d_mean_eucdist[d_mean_eucdist$tau.cat == "Low" | d_mean_eucdist$tau.cat == "Medium",]
d_mean_eucdist_lowhigh <- d_mean_eucdist[d_mean_eucdist$tau.cat == "Low" | d_mean_eucdist$tau.cat == "High",]
d_mean_eucdist_medhigh <- d_mean_eucdist[d_mean_eucdist$tau.cat == "Medium" | d_mean_eucdist$tau.cat == "High",]



plot_disttime_d.ls.lh <- ggplot(d_mean_eucdist_lowhigh, aes(x = gen, y = dist_mean, col = delmu.cat)) +
  facet_grid(locisigma.cat~tau.cat) +
  geom_line(position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Euclidean distance from optimum", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = t_lab
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_disttime_d.ls.lh)

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
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab_disttime_d.ls.lh <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab_disttime_d.ls.lh)




plot_disttime_d.ls.mh <- ggplot(d_mean_eucdist_medhigh, aes(x = gen, y = dist_mean, col = delmu.cat)) +
  facet_grid(locisigma.cat~tau.cat) +
  geom_line(position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Euclidean distance from optimum", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = t_lab
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_disttime_d.ls.mh)

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
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab_disttime_d.ls.mh <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab_disttime_d.ls.mh)




plot_disttime_d.ls.lm <- ggplot(d_mean_eucdist_lowmed, aes(x = gen, y = dist_mean, col = delmu.cat)) +
  facet_grid(locisigma.cat~tau.cat) +
  geom_line(position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Euclidean distance from optimum", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = t_lab
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_disttime_d.ls.lm)

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
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab_disttime_d.ls.lm <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab_disttime_d.ls.lm)

###########################################################
# Plot them together

library(gridExtra)


grid.arrange(plot_gtab_disttime.d.ls.n, plot_gtab_disttime_d.ls.lh, plot_gtab_disttime_d.ls.mh, plot_gtab_disttime_d.ls.lm, newpage = T)



# Altogether, much the same to read



plot_disttime_d.ls.s <- ggplot(d_mean_eucdist, aes(x = gen, y = dist_mean, col = delmu.cat)) +
  facet_grid(locisigma.cat~tau.cat) +
  geom_line(position = position_dodge(0.9), size = 0.1) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Euclidean distance from optimum", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = t_lab
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_disttime_d.ls.s)

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
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab_disttime_d.ls.s <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab_disttime_d.ls.s)

ggsave(filename = "disttime_d.ls.s.png", plot = plot_gtab_disttime_d.ls.s, width = 12, height = 8, dpi = 800)

###########################################################################
###########################################################################
# Simple figures: main effects of each parameter

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Delmu

dplot_eucdist_d <- d_eucdist_fingen[,c(4, 14)] %>%
  group_by(delmu.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_eucdist_d <- ggplot(dplot_eucdist_d, aes(x = delmu.cat, y = dist_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 500)) +
  theme_classic() +
  ggtitle(d_lab) +
  labs(x = d_lab, y = "\u03B4\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_eucdist_d

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# rwide

dplot_eucdist_r <- d_eucdist_fingen[,c(4, 16)] %>%
  group_by(rwide.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_eucdist_r <- ggplot(dplot_eucdist_r, aes(x = rwide.cat, y = dist_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 500)) +
  theme_classic() +
  ggtitle(r_lab) +
  labs(x = r_lab, y = "\u03B4\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_eucdist_r


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiorate

dplot_eucdist_pr <- d_eucdist_fingen[,c(4, 11)] %>%
  group_by(pleiorate.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_eucdist_pr <- ggplot(dplot_eucdist_pr, aes(x = pleiorate.cat, y = dist_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 500)) +
  theme_classic() +
  ggtitle(pr_lab) +
  labs(x = pr_lab, y = "\u03B4\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_eucdist_pr


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# locisigma

dplot_eucdist_ls <- d_eucdist_fingen[,c(4, 15)] %>%
  group_by(locisigma.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_eucdist_ls <- ggplot(dplot_eucdist_ls, aes(x = locisigma.cat, y = dist_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 500)) +
  theme_classic() +
  ggtitle(ls_lab) +
  labs(x = ls_lab, y = "\u03B4\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_eucdist_ls


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# tau

dplot_eucdist_t <- d_eucdist_fingen[,c(4, 13)] %>%
  group_by(tau.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_eucdist_t <- ggplot(dplot_eucdist_t, aes(x = tau.cat, y = dist_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 500)) +
  theme_classic() +
  ggtitle(t_lab) +
  labs(x = t_lab, y = "\u03B4\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        )

plot_eucdist_t

#############################################################################################################
# Graphs together
library(patchwork)

plot_dist_mainfx <- (plot_eucdist_d | plot_eucdist_r) / (plot_eucdist_ls | plot_eucdist_pr)
ggsave("plot_dist_mainfx.png", plot_dist_mainfx, height = )

plot_eucdist_t

# less ink: Title rather than axis 
